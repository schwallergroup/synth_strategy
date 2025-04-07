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
    This function detects a strategy involving the use of persistent directing groups
    (like methoxy) that remain throughout the synthesis.
    """
    # Track molecules with directing groups at each depth
    molecules_by_depth = {}
    max_depth = 0

    def dfs_traverse(node, depth=0, path_id=0):
        nonlocal max_depth
        max_depth = max(max_depth, depth)

        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            if mol_smiles:
                # Check specifically for methoxy groups (common directing groups)
                has_methoxy = "OC" in mol_smiles and checker.check_fg("Ether", mol_smiles)

                # Get more specific directing groups
                if has_methoxy:
                    # Store molecule info with the directing group
                    if depth not in molecules_by_depth:
                        molecules_by_depth[depth] = []

                    # Create RDKit molecule to analyze the methoxy positions
                    try:
                        mol = Chem.MolFromSmiles(mol_smiles)
                        if mol:
                            # Get atom indices of methoxy groups
                            methoxy_indices = []
                            for atom in mol.GetAtoms():
                                if atom.GetSymbol() == "O" and atom.GetDegree() == 2:
                                    neighbors = [n.GetAtomicNum() for n in atom.GetNeighbors()]
                                    if 6 in neighbors:  # Carbon
                                        for n in atom.GetNeighbors():
                                            if (
                                                n.GetAtomicNum() == 6 and n.GetDegree() == 1
                                            ):  # Methyl carbon
                                                methoxy_indices.append(atom.GetIdx())

                            molecules_by_depth[depth].append(
                                {
                                    "smiles": mol_smiles,
                                    "path_id": path_id,
                                    "methoxy_indices": methoxy_indices,
                                }
                            )
                            print(
                                f"Detected methoxy directing group at depth {depth}: {mol_smiles}"
                            )
                    except Exception as e:
                        print(f"Error analyzing molecule: {e}")

        # Traverse children with unique path IDs for each branch
        for i, child in enumerate(node.get("children", [])):
            dfs_traverse(child, depth + 1, path_id * 10 + i)

    # Start traversal
    dfs_traverse(route)

    # Check if we have any directing groups
    if not molecules_by_depth:
        print("No directing groups detected")
        return False

    # Check for persistence along synthesis paths
    persistent_paths = set()

    # Extract all unique path IDs at the earliest stage (highest depth)
    earliest_depth = max(molecules_by_depth.keys()) if molecules_by_depth else 0
    earliest_paths = set()
    for mol_info in molecules_by_depth.get(earliest_depth, []):
        earliest_paths.add(mol_info["path_id"])

    # For each path starting from the earliest stage, check if directing groups persist
    for path_id in earliest_paths:
        current_path_id = path_id
        consecutive_depths = []
        methoxy_positions = []

        # Trace the path from earliest to latest stage
        for depth in sorted(molecules_by_depth.keys(), reverse=True):
            for mol_info in molecules_by_depth[depth]:
                # Check if this molecule is in our current path
                if (
                    mol_info["path_id"] == current_path_id
                    or current_path_id % (10 * (10 ** (max_depth - depth))) == mol_info["path_id"]
                ):
                    consecutive_depths.append(depth)
                    methoxy_positions.append(mol_info["methoxy_indices"])
                    # Update path_id for next iteration (going up the tree)
                    current_path_id = mol_info["path_id"] // 10
                    break

        # Calculate coverage and check if methoxy persists
        if consecutive_depths:
            # Check if we have directing groups in at least 75% of the depths
            depth_coverage = len(consecutive_depths) / (max_depth + 1)

            # Check if methoxy groups persist (present in consecutive molecules)
            methoxy_persists = len(consecutive_depths) >= 3  # At least 3 consecutive depths

            if depth_coverage >= 0.5 and methoxy_persists:
                persistent_paths.add(path_id)
                print(
                    f"Path {path_id} has persistent methoxy directing group with coverage {depth_coverage:.2f}"
                )

    # Check if we have a continuous sequence of depths with methoxy groups
    depth_sequence = sorted(molecules_by_depth.keys())
    longest_sequence = 0
    current_sequence = 1

    for i in range(1, len(depth_sequence)):
        if depth_sequence[i] == depth_sequence[i - 1] + 2:  # Account for reaction nodes
            current_sequence += 1
        else:
            longest_sequence = max(longest_sequence, current_sequence)
            current_sequence = 1

    longest_sequence = max(longest_sequence, current_sequence)

    # Strategy is present if we have a long enough sequence of depths with methoxy groups
    # or if at least one path has persistent directing group
    strategy_present = (longest_sequence >= 3) or (len(persistent_paths) > 0)

    print(f"Depths with directing groups: {sorted(molecules_by_depth.keys())}")
    print(f"Max depth: {max_depth}")
    print(f"Longest sequence of consecutive depths: {longest_sequence}")
    print(f"Persistent paths: {persistent_paths}")
    print(f"Strategy detected: {strategy_present}")

    return strategy_present
