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

root_data = "/home/andres/Documents/steerable_retro/data"

fg_args = {
    "file_path": f"{root_data}/patterns/functional_groups.json",
    "value_field": "pattern",
    "key_field": "name",
}
reaction_class_args = {
    "file_path": f"{root_data}/patterns/smirks.json",
    "value_field": "smirks",
    "key_field": "name",
}
ring_smiles_args = {
    "file_path": f"{root_data}/patterns/chemical_rings_smiles.json",
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
    Detects if the synthesis involves transformation of a benzylic bromide to a methyl group.
    """
    transformation_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal transformation_detected

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            try:
                # Check if this is a relevant reaction type
                is_relevant_reaction = (
                    checker.check_reaction("Dehalogenation", rsmi)
                    or checker.check_reaction("Aromatic dehalogenation", rsmi)
                    or checker.check_reaction("Wohl-Ziegler bromination benzyl primary", rsmi)
                    or checker.check_reaction("Hydrogenolysis of tertiary amines", rsmi)
                )

                # Get reactant and product molecules
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]
                product_mol = Chem.MolFromSmiles(product)

                if product_mol and any(r for r in reactant_mols if r):
                    # Forward direction: benzylic bromide in reactants -> methyl in product
                    benzylic_br_in_reactants = False
                    benzylic_br_reactant_idx = -1

                    for r_idx, r_mol in enumerate(reactant_mols):
                        if (
                            r_mol
                            and checker.check_fg("Primary halide", reactants[r_idx])
                            and "Br" in reactants[r_idx]
                        ):
                            # Check if it's specifically a benzylic bromide
                            benzylic_br_pattern = Chem.MolFromSmarts("c[CH2]Br")
                            if r_mol.HasSubstructMatch(benzylic_br_pattern):
                                benzylic_br_in_reactants = True
                                benzylic_br_reactant_idx = r_idx
                                print(f"Found benzylic bromide in reactant: {reactants[r_idx]}")
                                break

                    # Check for benzylic methyl in product
                    benzylic_me_in_product = False
                    if product_mol:
                        benzylic_me_pattern = Chem.MolFromSmarts("c[CH3]")
                        if product_mol.HasSubstructMatch(benzylic_me_pattern):
                            benzylic_me_in_product = True
                            print(f"Found benzylic methyl in product: {product}")

                    # Reverse direction: benzylic methyl in reactants -> bromide in product (retrosynthetic)
                    benzylic_me_in_reactants = False
                    benzylic_me_reactant_idx = -1

                    for r_idx, r_mol in enumerate(reactant_mols):
                        if r_mol:
                            benzylic_me_pattern = Chem.MolFromSmarts("c[CH3]")
                            if r_mol.HasSubstructMatch(benzylic_me_pattern):
                                benzylic_me_in_reactants = True
                                benzylic_me_reactant_idx = r_idx
                                print(f"Found benzylic methyl in reactant: {reactants[r_idx]}")
                                break

                    # Check for benzylic bromide in product
                    benzylic_br_in_product = False
                    if (
                        product_mol
                        and checker.check_fg("Primary halide", product)
                        and "Br" in product
                    ):
                        benzylic_br_pattern = Chem.MolFromSmarts("c[CH2]Br")
                        if product_mol.HasSubstructMatch(benzylic_br_pattern):
                            benzylic_br_in_product = True
                            print(f"Found benzylic bromide in product: {product}")

                    # Verify transformation in either direction
                    if (benzylic_br_in_reactants and benzylic_me_in_product) or (
                        benzylic_me_in_reactants and benzylic_br_in_product
                    ):
                        # For atom-mapped reactions, try to ensure the same carbon is involved
                        has_matching_atom_mapping = False

                        # Check for any mapped carbon atom
                        for i in range(1, 100):  # Check reasonable range of atom mappings
                            atom_map = f"[C:{i}]"
                            if atom_map in rsmi:
                                if (
                                    benzylic_br_in_reactants
                                    and benzylic_me_in_product
                                    and atom_map in product
                                    and atom_map in reactants[benzylic_br_reactant_idx]
                                ):
                                    has_matching_atom_mapping = True
                                    break
                                elif (
                                    benzylic_me_in_reactants
                                    and benzylic_br_in_product
                                    and atom_map in product
                                    and atom_map in reactants[benzylic_me_reactant_idx]
                                ):
                                    has_matching_atom_mapping = True
                                    break

                        if has_matching_atom_mapping:
                            print(f"Confirmed transformation with atom mapping: {rsmi}")
                            transformation_detected = True
                        elif is_relevant_reaction:
                            print(f"Confirmed transformation with reaction type: {rsmi}")
                            transformation_detected = True
                        else:
                            # Additional check: if we can see the bromine is removed/added
                            if benzylic_br_in_reactants and benzylic_me_in_product:
                                r_mol = reactant_mols[benzylic_br_reactant_idx]
                                r_br_count = sum(
                                    1 for atom in r_mol.GetAtoms() if atom.GetSymbol() == "Br"
                                )
                                p_br_count = sum(
                                    1 for atom in product_mol.GetAtoms() if atom.GetSymbol() == "Br"
                                )

                                if r_br_count > p_br_count:
                                    print(
                                        f"Confirmed forward transformation by bromine count change: {rsmi}"
                                    )
                                    transformation_detected = True

                            elif benzylic_me_in_reactants and benzylic_br_in_product:
                                r_mol = reactant_mols[benzylic_me_reactant_idx]
                                r_br_count = sum(
                                    1 for atom in r_mol.GetAtoms() if atom.GetSymbol() == "Br"
                                )
                                p_br_count = sum(
                                    1 for atom in product_mol.GetAtoms() if atom.GetSymbol() == "Br"
                                )

                                if p_br_count > r_br_count:
                                    print(
                                        f"Confirmed reverse transformation by bromine count change: {rsmi}"
                                    )
                                    transformation_detected = True

                            # Check if this is a Wohl-Ziegler reaction in reverse
                            if checker.check_reaction(
                                "Wohl-Ziegler bromination benzyl primary", rsmi
                            ):
                                print(
                                    f"Confirmed transformation as Wohl-Ziegler bromination: {rsmi}"
                                )
                                transformation_detected = True
            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Benzylic bromide to methyl transformation: {transformation_detected}")
    return transformation_detected
