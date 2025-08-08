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
    This function detects a synthetic strategy involving a specific sequence of
    functional group interconversions: carboxylic acid → nitrile → halomethyl →
    hydroxymethyl → ester → carboxylic acid.
    """
    # Initialize tracking variables
    sequence = []
    visited_reactions = set()

    # Define the expected sequence
    expected_sequence = [
        "carboxylic_acid",
        "nitrile",
        "halomethyl",
        "hydroxymethyl",
        "ester",
        "carboxylic_acid",
    ]

    def dfs_traverse(node, depth=0):
        if (
            node["type"] == "reaction"
            and node.get("metadata", {}).get("rsmi")
            and node["metadata"]["rsmi"] not in visited_reactions
        ):
            rsmi = node["metadata"]["rsmi"]
            visited_reactions.add(rsmi)

            # Extract reactants and product
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for functional group transformations
            # Note: These are in forward reaction direction
            if checker.check_fg("Carboxylic acid", product) and any(
                checker.check_fg("Nitrile", r) for r in reactants
            ):
                # Nitrile to carboxylic acid
                sequence.append(("nitrile", "carboxylic_acid"))
                print(f"Found transformation: nitrile → carboxylic_acid at depth {depth}")

            elif checker.check_fg("Nitrile", product) and any(
                checker.check_fg("Primary halide", r) for r in reactants
            ):
                # Halomethyl to nitrile
                sequence.append(("halomethyl", "nitrile"))
                print(f"Found transformation: halomethyl → nitrile at depth {depth}")

            elif checker.check_fg("Primary halide", product) and any(
                checker.check_fg("Primary alcohol", r) for r in reactants
            ):
                # Hydroxymethyl to halomethyl (corrected direction)
                sequence.append(("hydroxymethyl", "halomethyl"))
                print(f"Found transformation: hydroxymethyl → halomethyl at depth {depth}")

            elif checker.check_fg("Primary alcohol", product) and any(
                checker.check_fg("Ester", r) for r in reactants
            ):
                # Ester to hydroxymethyl
                sequence.append(("ester", "hydroxymethyl"))
                print(f"Found transformation: ester → hydroxymethyl at depth {depth}")

            elif checker.check_fg("Ester", product) and any(
                checker.check_fg("Carboxylic acid", r) for r in reactants
            ):
                # Carboxylic acid to ester
                sequence.append(("carboxylic_acid", "ester"))
                print(f"Found transformation: carboxylic_acid → ester at depth {depth}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Transformations found: {sequence}")

    # Build a directed graph of transformations
    graph = {}
    for src, dst in sequence:
        if src not in graph:
            graph[src] = []
        graph[src].append(dst)

    # Find all possible paths in the sequence
    def find_paths(graph, start, path=None):
        if path is None:
            path = [start]

        paths = [path]

        if start in graph:
            for next_node in graph[start]:
                if next_node not in path:  # Avoid cycles
                    new_paths = find_paths(graph, next_node, path + [next_node])
                    paths.extend(new_paths)

        return paths

    # Find all paths starting from each node
    all_paths = []
    for node in set([src for src, _ in sequence] + [dst for _, dst in sequence]):
        paths = find_paths(graph, node)
        all_paths.extend(paths)

    # Filter paths to only include those with at least 3 nodes
    valid_paths = [p for p in all_paths if len(p) >= 3]

    # Check if any valid path is a subsequence of the expected sequence
    for path in valid_paths:
        # Check if path is a subsequence of expected_sequence
        i, j = 0, 0
        while i < len(expected_sequence) and j < len(path):
            if expected_sequence[i] == path[j]:
                j += 1
            i += 1

        if j == len(path):  # All elements of path were found in expected_sequence in order
            print(f"Found valid path: {path}")

            # Check if the path contains at least 3 consecutive elements from expected_sequence
            for i in range(len(expected_sequence) - 2):
                subseq = expected_sequence[i : i + 3]

                # Check if subseq is a subsequence of path
                sub_i, sub_j = 0, 0
                while sub_i < len(path) and sub_j < len(subseq):
                    if path[sub_i] == subseq[sub_j]:
                        sub_j += 1
                    sub_i += 1

                if sub_j == len(subseq):  # All elements of subseq were found in path in order
                    print(f"Found valid subsequence: {subseq}")
                    return True

    # Check if we have at least 3 consecutive steps from the expected sequence
    # Convert sequence to a list of functional groups
    fg_sequence = []
    for src, dst in sequence:
        if not fg_sequence:
            fg_sequence.append(src)
        fg_sequence.append(dst)

    # Remove duplicates while preserving order
    seen = set()
    unique_fg_sequence = []
    for fg in fg_sequence:
        if fg not in seen:
            seen.add(fg)
            unique_fg_sequence.append(fg)

    print(f"Functional group sequence: {unique_fg_sequence}")

    # Check for consecutive subsequences of length 3 or more
    if len(unique_fg_sequence) >= 3:
        for i in range(len(expected_sequence) - 2):
            window = expected_sequence[i : i + 3]

            # Check if window is a subsequence of unique_fg_sequence
            sub_i, sub_j = 0, 0
            while sub_i < len(unique_fg_sequence) and sub_j < len(window):
                if unique_fg_sequence[sub_i] == window[sub_j]:
                    sub_j += 1
                sub_i += 1

            if sub_j == len(
                window
            ):  # All elements of window were found in unique_fg_sequence in order
                print(f"Found valid subsequence: {window}")
                return True

    return False
