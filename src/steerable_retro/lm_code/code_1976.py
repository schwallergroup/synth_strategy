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
    Detects if the synthesis route has a late-stage esterification
    (ester formation in the final step or penultimate steps)
    """
    late_stage_has_esterification = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_has_esterification

        # Consider final step and penultimate steps as "late-stage"
        if node["type"] == "reaction" and depth <= 2:  # Final or near-final steps
            print(f"Examining reaction step at depth {depth}")
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                print(f"Reaction SMILES: {rsmi}")

                # Check for various esterification reactions
                esterification_reactions = [
                    "Esterification of Carboxylic Acids",
                    "Schotten-Baumann to ester",
                    "Transesterification",
                    "O-alkylation of carboxylic acids with diazo compounds",
                    "Oxidative esterification of primary alcohols",
                ]

                for rxn_type in esterification_reactions:
                    if checker.check_reaction(rxn_type, rsmi):
                        print(f"Detected {rxn_type} reaction at depth {depth}")
                        late_stage_has_esterification = True
                        return

                # Alternative check: look for ester formation
                try:
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    # Create molecule objects
                    product_mol = Chem.MolFromSmiles(product)
                    reactant_mols = [
                        Chem.MolFromSmiles(r) for r in reactants if Chem.MolFromSmiles(r)
                    ]

                    if product_mol:
                        # Check for ester group in product
                        ester_pattern = Chem.MolFromSmarts("C(=O)O[#6]")
                        carbonate_pattern = Chem.MolFromSmarts("O-C(=O)-O-[#6]")

                        product_has_ester = product_mol.HasSubstructMatch(ester_pattern)
                        product_has_carbonate = product_mol.HasSubstructMatch(carbonate_pattern)

                        if product_has_ester or product_has_carbonate:
                            print(
                                f"Product contains {'carbonate' if product_has_carbonate else 'ester'} group"
                            )

                            # Check if ester/carbonate was formed in this step (not present in reactants)
                            ester_in_reactants = any(
                                mol.HasSubstructMatch(ester_pattern) for mol in reactant_mols if mol
                            )
                            carbonate_in_reactants = any(
                                mol.HasSubstructMatch(carbonate_pattern)
                                for mol in reactant_mols
                                if mol
                            )

                            # Check for esterification reagents
                            acid_pattern = Chem.MolFromSmarts("C(=O)O[H]")
                            acyl_halide_pattern = Chem.MolFromSmarts("C(=O)[F,Cl,Br,I]")
                            chloroformate_pattern = Chem.MolFromSmarts("[F,Cl,Br,I]C(=O)O[#6]")
                            anhydride_pattern = Chem.MolFromSmarts("C(=O)OC(=O)")

                            acid_in_reactants = any(
                                mol.HasSubstructMatch(acid_pattern) for mol in reactant_mols if mol
                            )
                            acyl_halide_in_reactants = any(
                                mol.HasSubstructMatch(acyl_halide_pattern)
                                for mol in reactant_mols
                                if mol
                            )
                            chloroformate_in_reactants = any(
                                mol.HasSubstructMatch(chloroformate_pattern)
                                for mol in reactant_mols
                                if mol
                            )
                            anhydride_in_reactants = any(
                                mol.HasSubstructMatch(anhydride_pattern)
                                for mol in reactant_mols
                                if mol
                            )

                            # Check for alcohol/phenol in reactants
                            alcohol_pattern = Chem.MolFromSmarts("[#6]-[O;H1]")
                            phenol_pattern = Chem.MolFromSmarts("c-[O;H1]")

                            alcohol_in_reactants = any(
                                mol.HasSubstructMatch(alcohol_pattern)
                                for mol in reactant_mols
                                if mol
                            )
                            phenol_in_reactants = any(
                                mol.HasSubstructMatch(phenol_pattern)
                                for mol in reactant_mols
                                if mol
                            )

                            # Detect esterification
                            if (
                                (not ester_in_reactants and not carbonate_in_reactants)
                                and (
                                    (
                                        acid_in_reactants
                                        or acyl_halide_in_reactants
                                        or anhydride_in_reactants
                                    )
                                    and (alcohol_in_reactants or phenol_in_reactants)
                                )
                                or (
                                    chloroformate_in_reactants
                                    and (alcohol_in_reactants or phenol_in_reactants)
                                )
                            ):
                                print(f"Ester/carbonate was formed in step at depth {depth}")
                                late_stage_has_esterification = True
                            elif ester_in_reactants or carbonate_in_reactants:
                                # Check for transesterification
                                if checker.check_reaction("Transesterification", rsmi):
                                    print(f"Detected transesterification at depth {depth}")
                                    late_stage_has_esterification = True
                        else:
                            print("Product does not contain ester or carbonate group")
                except Exception as e:
                    print(f"Error processing reaction SMILES: {e}")
            else:
                print("No reaction SMILES found in metadata")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    print(f"Final result: {late_stage_has_esterification}")
    return late_stage_has_esterification
