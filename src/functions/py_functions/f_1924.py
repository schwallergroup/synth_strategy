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
    This function detects if the synthesis follows a convergent approach with 3 or more fragments
    coming together in the final steps.
    """
    # Track all complex reactants at each depth and by branch
    reactants_by_depth = {}
    branches = {}
    branch_id = 0

    def dfs_traverse(node, depth=0, branch=0):
        nonlocal branch_id

        if node["type"] == "reaction":
            try:
                # Extract reactants
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")

                # Identify complex reactants (with >8 heavy atoms or containing rings)
                complex_reactants = []
                for r in reactants_smiles:
                    mol = Chem.MolFromSmiles(r)
                    if mol:
                        # Consider a reactant complex if it has >8 heavy atoms or contains a ring
                        is_complex = mol.GetNumHeavyAtoms() > 8

                        # Check for presence of rings
                        ring_info = mol.GetRingInfo()
                        has_ring = ring_info.NumRings() > 0

                        # Check for presence of functional groups that indicate complexity
                        has_complex_fg = False
                        for fg in [
                            "Aromatic halide",
                            "Ester",
                            "Amide",
                            "Nitrile",
                            "Carboxylic acid",
                        ]:
                            if checker.check_fg(fg, r):
                                has_complex_fg = True
                                break

                        if is_complex or has_ring or has_complex_fg:
                            complex_reactants.append(r)

                # Store all complex reactants at this depth
                if complex_reactants:
                    reactants_by_depth.setdefault(depth, []).extend(complex_reactants)

                    # If we have multiple complex reactants, create new branches for them
                    if len(complex_reactants) > 1:
                        for i, r in enumerate(complex_reactants):
                            if i > 0:  # First reactant continues current branch
                                branch_id += 1
                                branches.setdefault(branch_id, []).append(r)
                            else:
                                branches.setdefault(branch, []).append(r)
                    else:
                        branches.setdefault(branch, []).extend(complex_reactants)

                print(
                    f"Depth {depth}: Found {len(complex_reactants)} complex reactants out of {len(reactants_smiles)} total"
                )
                for i, r in enumerate(reactants_smiles):
                    mol = Chem.MolFromSmiles(r)
                    if mol:
                        heavy_atoms = mol.GetNumHeavyAtoms()
                        print(f"  Reactant {i+1}: {r} - {heavy_atoms} heavy atoms")
            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        # Continue traversing with multiple branches if needed
        if len(node.get("children", [])) > 1:
            # Create new branches for multiple children
            for i, child in enumerate(node.get("children", [])):
                if i > 0:  # First child continues current branch
                    branch_id += 1
                    dfs_traverse(child, depth + 1, branch_id)
                else:
                    dfs_traverse(child, depth + 1, branch)
        else:
            # Continue with same branch
            for child in node.get("children", []):
                dfs_traverse(child, depth + 1, branch)

    # Start traversal
    dfs_traverse(route)

    # Check if any of the final 4 steps (depths 0-3) have 3 or more unique complex reactants
    max_depth = max(reactants_by_depth.keys()) if reactants_by_depth else 0
    check_depths = range(min(4, max_depth + 1))

    for depth in check_depths:
        if depth in reactants_by_depth:
            # Count unique complex reactants at this depth
            unique_reactants = set(reactants_by_depth[depth])
            if len(unique_reactants) >= 3:
                print(
                    f"Detected convergent synthesis with {len(unique_reactants)} unique complex fragments at depth {depth}"
                )
                print(f"Complex fragments: {list(unique_reactants)}")
                return True

    # Check if depths 0-3 combined have 3 or more unique complex reactants
    all_late_stage_reactants = set()
    for depth in check_depths:
        if depth in reactants_by_depth:
            all_late_stage_reactants.update(reactants_by_depth[depth])

    if len(all_late_stage_reactants) >= 3:
        print(
            f"Detected convergent synthesis with {len(all_late_stage_reactants)} unique complex fragments across the final {len(check_depths)} steps"
        )
        print(f"Complex fragments: {list(all_late_stage_reactants)}")
        return True

    # Check if we have at least 3 distinct branches with complex reactants
    if len(branches) >= 3:
        print(f"Detected convergent synthesis with {len(branches)} distinct branches")
        return True

    # Count total unique complex reactants across all depths
    all_reactants = set()
    for reactants in reactants_by_depth.values():
        all_reactants.update(reactants)

    if len(all_reactants) >= 3:
        print(
            f"Detected convergent synthesis with {len(all_reactants)} total unique complex fragments"
        )
        return True

    print("No convergent synthesis pattern detected")
    return False
