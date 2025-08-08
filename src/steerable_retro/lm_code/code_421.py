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
    Detects if the morpholine ring is preserved throughout the synthesis.
    This means that once the morpholine ring is introduced, it remains intact
    in all subsequent steps leading to the final product.
    """
    # Track when morpholine is first introduced and if it's preserved
    morpholine_introduced = False
    morpholine_preserved = True

    def dfs_traverse(node, depth=0, parent_has_morpholine=False):
        nonlocal morpholine_introduced, morpholine_preserved

        if node["type"] == "mol" and "smiles" in node:
            # Check if this molecule has a morpholine ring
            has_morpholine = checker.check_ring("morpholine", node["smiles"])

            # Skip checking in-stock molecules (starting materials)
            if node.get("in_stock", False):
                print(f"Skipping in-stock molecule at depth {depth}: {node['smiles']}")
                return

            # If parent had morpholine but this molecule doesn't, preservation failed
            if parent_has_morpholine and not has_morpholine:
                morpholine_preserved = False
                print(f"Morpholine lost at depth {depth} in molecule {node['smiles']}")

            # If this molecule has morpholine, mark it as introduced
            if has_morpholine:
                morpholine_introduced = True
                print(f"Morpholine found at depth {depth} in molecule {node['smiles']}")

        # Continue traversing with updated parent_has_morpholine status
        current_has_morpholine = False
        if node["type"] == "mol" and "smiles" in node:
            current_has_morpholine = checker.check_ring("morpholine", node["smiles"])

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1, current_has_morpholine)

    # First check if the target molecule contains morpholine
    target_has_morpholine = checker.check_ring("morpholine", route["smiles"])
    if not target_has_morpholine:
        print(f"Target molecule does not contain morpholine: {route['smiles']}")
        return False

    # Start traversal from the root
    dfs_traverse(route)

    # Return True if morpholine was introduced and preserved throughout
    if morpholine_introduced and morpholine_preserved:
        print("Morpholine preserved throughout synthesis after introduction")
        return True
    elif not morpholine_introduced:
        print("Morpholine not found in any non-stock molecules in the route")
        return False
    else:
        print("Morpholine was introduced but not preserved throughout synthesis")
        return False
