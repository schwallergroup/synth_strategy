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
    This function detects if the synthesis route includes a late-stage amide to nitrile conversion.
    Late stage means in the first half of the synthesis (low depth in retrosynthetic tree).
    """
    # First pass: calculate the maximum depth
    max_depth = 0

    def calculate_max_depth(node, depth=0):
        nonlocal max_depth
        max_depth = max(max_depth, depth)
        for child in node.get("children", []):
            calculate_max_depth(child, depth + 1)

    calculate_max_depth(route)
    print(f"Maximum depth of synthesis route: {max_depth}")

    # Second pass: detect amide to nitrile conversion
    amide_to_nitrile_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal amide_to_nitrile_detected

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if any reactant contains amide and product contains nitrile
            has_amide = False
            for reactant in reactants:
                if (
                    checker.check_fg("Primary amide", reactant)
                    or checker.check_fg("Secondary amide", reactant)
                    or checker.check_fg("Tertiary amide", reactant)
                ):
                    has_amide = True
                    print(f"Found amide in reactant: {reactant}")
                    break

            has_nitrile = checker.check_fg("Nitrile", product)
            if has_nitrile:
                print(f"Found nitrile in product: {product}")

            # Check if this is an amide to nitrile conversion
            if has_amide and has_nitrile:
                print(f"Potential amide to nitrile conversion at depth {depth}")

                # Check for dehydration reactions that could convert amide to nitrile
                is_dehydration = False

                # Try to check if this is a dehydration reaction (amide to nitrile)
                # This is typically done with dehydrating agents like POCl3, SOCl2, P2O5, etc.
                if checker.check_reaction("Dehydration", rsmi):
                    is_dehydration = True
                    print("Confirmed dehydration reaction")

                # If we can't identify the specific reaction type, check if the reaction
                # involves removal of water/oxygen which is characteristic of amide to nitrile conversion
                if not is_dehydration:
                    # Count O atoms in reactants and products to check for oxygen loss
                    reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]
                    product_mol = Chem.MolFromSmiles(product)

                    if all(reactant_mols) and product_mol:
                        reactant_o_count = sum(
                            atom.GetSymbol() == "O"
                            for mol in reactant_mols
                            for atom in mol.GetAtoms()
                        )
                        product_o_count = sum(
                            atom.GetSymbol() == "O" for atom in product_mol.GetAtoms()
                        )

                        if reactant_o_count > product_o_count:
                            is_dehydration = True
                            print(
                                f"Oxygen loss detected: {reactant_o_count} in reactants, {product_o_count} in product"
                            )

                if is_dehydration:
                    # Only consider it late-stage if in first half of synthesis
                    if depth < max_depth / 2:
                        print(
                            f"Late-stage amide to nitrile conversion detected at depth {depth} (max depth: {max_depth})"
                        )
                        amide_to_nitrile_detected = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    print(f"Final result: {amide_to_nitrile_detected}")

    return amide_to_nitrile_detected
