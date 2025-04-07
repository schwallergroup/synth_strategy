#!/bin/python

"""LM-defined function for strategy description."""

import copy
import re
from collections import deque

import rdkit
import rdkit.Chem as Chem
from rdkit import Chem
from rdkit.Chem import (
    AllChem,
    Descriptors,
    Lipinski,
    rdChemReactions,
    rdFMCS,
    rdMolDescriptors,
    rdmolops,
)
from rdkit.Chem.Scaffolds import MurckoScaffold
from steerable_retro.utils import check, fuzzy_dict
from steerable_retro.utils.check import Check

fg_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/functional_groups.json",
    "value_field": "pattern",
    "key_field": "name",
}
reaction_class_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/smirks.json",
    "value_field": "smirks",
    "key_field": "name",
}
ring_smiles_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/chemical_rings_smiles.json",
    "value_field": "smiles",
    "key_field": "name",
}
functional_groups = fuzzy_dict.FuzzyDict.from_json(**fg_args)
reaction_classes = fuzzy_dict.FuzzyDict.from_json(**reaction_class_args)
ring_smiles = fuzzy_dict.FuzzyDict.from_json(**ring_smiles_args)

checker = check.Check(
    fg_dict=functional_groups, reaction_dict=reaction_classes, ring_dict=ring_smiles
)


def main(route):
    """
    This function detects a synthetic strategy involving multiple sequential N-methylation
    reactions on a nitrogen-rich heterocycle, specifically with an imidazolidinedione scaffold.
    """
    # Track the number of N-methylation reactions
    n_methylation_count = 0
    # Flag to check if imidazolidinedione scaffold is present
    has_imidazolidinedione = False
    # Track molecules that have undergone N-methylation
    methylated_molecules = {}  # Store molecule SMILES and their methylation status

    def is_n_methylation(rxn_smiles):
        """Check if a reaction is an N-methylation using the checker function"""
        try:
            # Check for various methylation reaction types
            return (
                checker.check_reaction("N-methylation", rxn_smiles)
                or checker.check_reaction("Methylation", rxn_smiles)
                or checker.check_reaction("Methylation with MeI_primary", rxn_smiles)
                or checker.check_reaction("Methylation with MeI_secondary", rxn_smiles)
                or checker.check_reaction("Methylation with MeI_tertiary", rxn_smiles)
                or checker.check_reaction("DMS Amine methylation", rxn_smiles)
                or checker.check_reaction(
                    "Eschweiler-Clarke Primary Amine Methylation", rxn_smiles
                )
                or checker.check_reaction(
                    "Eschweiler-Clarke Secondary Amine Methylation", rxn_smiles
                )
                or checker.check_reaction(
                    "Reductive methylation of primary amine with formaldehyde",
                    rxn_smiles,
                )
            )
        except Exception as e:
            print(f"Error in is_n_methylation: {e}")
            return False

    def check_imidazolidinedione(smiles):
        """Check if a molecule contains an imidazolidinedione scaffold (hydantoin)"""
        try:
            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                return False

            # Check using the ring pattern directly
            if checker.check_ring("imidazolidine", smiles):
                # Now check if it has the dione pattern
                pattern = Chem.MolFromSmarts("O=CN[C,N]C(=O)N")
                if mol.HasSubstructMatch(pattern):
                    return True

            # SMARTS pattern for imidazolidinedione (hydantoin)
            pattern = Chem.MolFromSmarts("O=C1NC(=O)[C,N]N1")
            if mol.HasSubstructMatch(pattern):
                return True

            # Alternative pattern for N-methylated hydantoin
            pattern2 = Chem.MolFromSmarts("O=C1N(C)C(=O)[C,N]N1")
            if mol.HasSubstructMatch(pattern2):
                return True

            # Alternative pattern for di-N-methylated hydantoin
            pattern3 = Chem.MolFromSmarts("O=C1N(C)C(=O)[C,N]N1C")
            if mol.HasSubstructMatch(pattern3):
                return True

            return False
        except Exception as e:
            print(f"Error in check_imidazolidinedione: {e}")
            return False

    def count_n_methyl_groups(smiles):
        """Count the number of N-methyl groups in the imidazolidinedione scaffold"""
        try:
            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                return 0

            # Pattern for N-methyl in hydantoin at position 1
            pattern1 = Chem.MolFromSmarts("CN1C(=O)NC(=O)[C,N]1")
            # Pattern for N-methyl in hydantoin at position 3
            pattern2 = Chem.MolFromSmarts("C1NC(=O)N(C)C(=O)[C,N]1")
            # Pattern for di-N-methylated hydantoin
            pattern3 = Chem.MolFromSmarts("CN1C(=O)N(C)C(=O)[C,N]1")

            count = 0
            if mol.HasSubstructMatch(pattern1):
                count += 1
            if mol.HasSubstructMatch(pattern2):
                count += 1
            # If we have a di-methylated pattern, ensure we count it as 2
            if mol.HasSubstructMatch(pattern3):
                count = 2

            return count
        except Exception as e:
            print(f"Error in count_n_methyl_groups: {e}")
            return 0

    def is_methylation_change(product_smiles, reactant_smiles):
        """Check if there's an increase in N-methyl groups from reactant to product"""
        try:
            product_methyl_count = count_n_methyl_groups(product_smiles)
            reactant_methyl_count = count_n_methyl_groups(reactant_smiles)

            print(
                f"Methylation change check: {reactant_smiles} ({reactant_methyl_count} methyl groups) -> {product_smiles} ({product_methyl_count} methyl groups)"
            )
            return product_methyl_count > reactant_methyl_count
        except Exception as e:
            print(f"Error in is_methylation_change: {e}")
            return False

    def dfs_traverse(node, depth=0, parent_smiles=None):
        nonlocal n_methylation_count, has_imidazolidinedione

        if node["type"] == "mol":
            mol_smiles = node["smiles"]

            # Check if molecule contains imidazolidinedione scaffold
            if check_imidazolidinedione(mol_smiles):
                has_imidazolidinedione = True
                print(f"Found imidazolidinedione scaffold in molecule: {mol_smiles}")

                # Store the methylation status of this molecule
                methyl_count = count_n_methyl_groups(mol_smiles)
                methylated_molecules[mol_smiles] = methyl_count
                print(
                    f"Molecule has {methyl_count} N-methyl groups on imidazolidinedione"
                )

        elif node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if this is an N-methylation reaction
                if is_n_methylation(rsmi):
                    print(f"Found potential N-methylation reaction: {rsmi}")

                    # Check if the product contains imidazolidinedione scaffold
                    if check_imidazolidinedione(product):
                        # Check if any reactant also contains imidazolidinedione
                        for reactant in reactants:
                            if check_imidazolidinedione(reactant):
                                # Verify that N-methylation occurred on the imidazolidinedione
                                if is_methylation_change(product, reactant):
                                    n_methylation_count += 1
                                    print(
                                        f"Confirmed N-methylation on imidazolidinedione: {reactant} -> {product}"
                                    )

                                    # Store the product's methylation state for future reference
                                    methylated_molecules[
                                        product
                                    ] = count_n_methyl_groups(product)

                                    # Check for sequential methylation
                                    if reactant in methylated_molecules:
                                        print(
                                            f"Sequential methylation detected: {reactant} ({methylated_molecules[reactant]} methyl groups) -> {product} ({count_n_methyl_groups(product)} methyl groups)"
                                        )
            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1, node.get("smiles", ""))

    # Start traversal
    dfs_traverse(route)

    # Strategy criteria: at least 2 N-methylation reactions and presence of imidazolidinedione scaffold
    strategy_detected = n_methylation_count >= 2 and has_imidazolidinedione
    print(f"Sequential N-methylation strategy detected: {strategy_detected}")
    print(f"N-methylation count: {n_methylation_count}")
    print(f"Imidazolidinedione scaffold present: {has_imidazolidinedione}")

    # Alternative detection method: check if we have molecules with different methylation states
    if has_imidazolidinedione and len(methylated_molecules) >= 2:
        methyl_counts = list(methylated_molecules.values())
        if 0 in methyl_counts and (1 in methyl_counts or 2 in methyl_counts):
            print("Detected sequential methylation based on molecule analysis")
            strategy_detected = True

    return strategy_detected
