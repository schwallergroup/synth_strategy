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
    Detects if an aldehyde functional group is maintained throughout the synthesis.
    This means the aldehyde group must be present in the final product and preserved
    through multiple reaction steps.
    """
    # Track paths where aldehyde is maintained
    aldehyde_paths = []

    def dfs_traverse(node, depth=0, path=None):
        if path is None:
            path = []

        current_path = path.copy()

        # For molecule nodes, check if aldehyde is present
        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            has_aldehyde = checker.check_fg("Aldehyde", mol_smiles)

            if has_aldehyde:
                current_path.append((depth, "mol", mol_smiles))
                print(f"Found aldehyde in molecule at depth {depth}: {mol_smiles}")
            elif current_path:  # If we had an aldehyde path but lost it
                # Save the path if it's long enough
                if len(current_path) >= 3:
                    aldehyde_paths.append(current_path)
                current_path = []  # Reset path since aldehyde is not maintained

        # For reaction nodes, check if aldehyde is maintained through the reaction
        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if aldehyde is in both reactants and product
            reactant_has_aldehyde = any(checker.check_fg("Aldehyde", r) for r in reactants)
            product_has_aldehyde = checker.check_fg("Aldehyde", product)

            if reactant_has_aldehyde and product_has_aldehyde:
                current_path.append((depth, "reaction", rsmi))
                print(f"Aldehyde maintained through reaction at depth {depth}: {rsmi}")
            elif current_path:  # If we had an aldehyde path but lost it
                # Save the path if it's long enough
                if len(current_path) >= 3:
                    aldehyde_paths.append(current_path)
                current_path = []  # Reset path since aldehyde is not maintained

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1, current_path)

        # If we're at the root and have a valid path, save it
        if depth == 0 and current_path and len(current_path) >= 3:
            aldehyde_paths.append(current_path)

    # Start traversal
    dfs_traverse(route)

    # Check if we found any valid paths where aldehyde is maintained
    if aldehyde_paths:
        longest_path = max(aldehyde_paths, key=len)
        depths = [d for d, _, _ in longest_path]
        print(f"Aldehyde maintained throughout synthesis at depths: {depths}")
        return True

    return False
