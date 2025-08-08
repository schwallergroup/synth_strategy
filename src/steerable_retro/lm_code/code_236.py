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
    This function detects if a nitrile group is preserved through multiple steps of the synthesis.
    """
    nitrile_by_depth = {}

    def dfs_traverse(node, current_depth=0):
        if node["type"] == "mol":
            smiles = node["smiles"]

            # Check for nitrile group using the checker function
            has_nitrile = checker.check_fg("Nitrile", smiles)
            print(
                f"Checking nitrile at depth {current_depth}, SMILES: {smiles}, Result: {has_nitrile}"
            )

            # Store nitrile presence at this depth
            if current_depth not in nitrile_by_depth:
                nitrile_by_depth[current_depth] = has_nitrile
            else:
                # If we already have a value for this depth, OR it with the new value
                # This handles branching synthesis routes
                nitrile_by_depth[current_depth] = nitrile_by_depth[current_depth] or has_nitrile

            if has_nitrile:
                print(f"Found nitrile group at depth {current_depth}")

        # Traverse children with incremented depth
        for child in node.get("children", []):
            dfs_traverse(child, current_depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Count how many depths have nitrile groups
    nitrile_depths = sorted([d for d in nitrile_by_depth.keys() if nitrile_by_depth[d]])

    # Check if nitrile is preserved through multiple steps (at least 2 different depths)
    nitrile_preservation = len(nitrile_depths) >= 2

    print(f"Nitrile preservation strategy detected: {nitrile_preservation}")
    print(f"Nitrile found at depths: {nitrile_depths}")

    return nitrile_preservation
