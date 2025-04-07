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
    This function detects if key functional groups (nitro-aromatic, allyl, dimethylamide)
    are maintained throughout the synthesis.
    """
    # Define the functional groups to track
    target_functional_groups = ["Nitro group", "Allyl", "Tertiary amide"]

    # Track which functional groups are present in the target molecule
    target_fg_present = {fg: False for fg in target_functional_groups}

    # Process the target molecule first (depth 0)
    if route["type"] == "mol":
        mol_smiles = route["smiles"]
        for fg in target_functional_groups:
            target_fg_present[fg] = checker.check_fg(fg, mol_smiles)
            print(f"Target molecule has {fg}: {target_fg_present[fg]}")

    # Store complete paths where each functional group persists
    complete_paths = []

    # Track functional groups through each path
    def dfs_traverse(node, current_path, path_fg_present, depth=0):
        # Make a copy of the current state to avoid modifying the parent's state
        current_path = current_path.copy()
        path_fg_present = copy.deepcopy(path_fg_present)

        # Add current node to path
        current_path.append(node)

        if node["type"] == "mol":
            mol_smiles = node["smiles"]

            # Check which functional groups are present in this molecule
            for fg in target_functional_groups:
                has_fg = checker.check_fg(fg, mol_smiles)

                # Update the current path state - only track FGs that were in target
                if target_fg_present[fg]:
                    path_fg_present[fg] = has_fg

            # If this is a starting material (leaf node), we've completed a path
            if node.get("in_stock", False) or not node.get("children", []):
                # Store this complete path and its FG persistence
                complete_paths.append((current_path, path_fg_present))
                print(f"Complete path found, FGs at depth {depth}: {path_fg_present}")
                return

        elif node["type"] == "reaction":
            # For reaction nodes, check if the reaction preserves our functional groups
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # In retrosynthesis, we're going from product to reactants
                # Check if FGs in reactants were present in the product
                for fg in target_functional_groups:
                    # Only check FGs that we're tracking (were in target)
                    if target_fg_present[fg]:
                        # Check if any reactant has this FG
                        fg_in_reactants = any(checker.check_fg(fg, r) for r in reactants_smiles)

                        # If FG is in reactants but not in product, it's created in this reaction
                        # (going backward), so not persistent from earlier
                        if fg_in_reactants and not checker.check_fg(fg, product_smiles):
                            print(f"FG {fg} created in reaction (not persistent): {rsmi}")
                            path_fg_present[fg] = False
            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Continue traversing children
        for child in node.get("children", []):
            dfs_traverse(child, current_path, path_fg_present, depth + 1)

    # Start traversal from the root
    dfs_traverse(route, [], target_fg_present)

    # Count functional groups that persist in at least one complete path
    persistent_fgs = set()
    for path, fg_status in complete_paths:
        for fg in target_functional_groups:
            if target_fg_present[fg] and fg_status[fg]:
                persistent_fgs.add(fg)

    persistent_fg_count = len(persistent_fgs)
    print(f"Target FGs present: {target_fg_present}")
    print(f"Persistent functional groups: {persistent_fgs}")
    print(f"Persistent functional groups count: {persistent_fg_count}")

    # Return True if at least two functional groups persist throughout
    return persistent_fg_count >= 2
