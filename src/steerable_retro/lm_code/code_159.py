#!/bin/python

"""LM-defined function for strategy description."""

from rdkit.Chem import AllChem, rdFMCS
import copy
from collections import deque
import rdkit.Chem as Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdChemReactions
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
import rdkit.Chem.rdFMCS
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem import AllChem, Descriptors, Lipinski
from rdkit.Chem import rdmolops
import re
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import AllChem, Descriptors
import traceback
import rdkit
from collections import Counter
from steerable_retro.utils.check import Check
from steerable_retro.utils import fuzzy_dict, check

root_data = "/home/dparm/steerable_retro/data"

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
    This function detects if an oxetane ring is present throughout the synthesis.
    The oxetane must be maintained in every step of the synthesis route.
    """
    # Track complete paths with oxetane
    paths_with_oxetane = []

    def dfs_traverse(node, current_path=None, depth=0):
        if current_path is None:
            current_path = []

        # Create a copy of the current path
        path = current_path.copy()

        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            has_oxetane = checker.check_ring("oxetane", mol_smiles)

            # Add molecule info to the path
            path.append(
                {"depth": depth, "smiles": mol_smiles, "has_oxetane": has_oxetane, "type": "mol"}
            )

            print(f"Depth {depth}: {'Found' if has_oxetane else 'No'} oxetane ring in {mol_smiles}")

            # If this is a leaf node (no children), save the complete path
            if not node.get("children", []):
                paths_with_oxetane.append(path)
        else:  # Reaction node
            # Add reaction info to the path
            reaction_smiles = node.get("metadata", {}).get("rsmi", "")
            path.append({"depth": depth, "rsmi": reaction_smiles, "type": "reaction"})

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, path, depth + 1)

    # Start traversal
    dfs_traverse(route)

    if not paths_with_oxetane:
        print("No complete paths found")
        return False

    # Check each complete path
    valid_paths = []
    for path in paths_with_oxetane:
        # Filter out reaction nodes for oxetane check
        mol_nodes = [node for node in path if node["type"] == "mol"]

        # Check if all molecule nodes in this path have oxetane
        if all(node["has_oxetane"] for node in mol_nodes):
            # Check if path has at least 3 steps (2 reactions)
            if len(mol_nodes) >= 3:
                valid_paths.append(path)
                print(f"Valid path found with {len(mol_nodes)} molecules, all containing oxetane")
            else:
                print(f"Path too short: only {len(mol_nodes)} molecules")
        else:
            print(f"Path invalid: not all molecules contain oxetane")

    # Return True if at least one valid path exists
    result = len(valid_paths) > 0
    print(f"Found {len(valid_paths)} valid paths with oxetane maintained throughout")
    return result
