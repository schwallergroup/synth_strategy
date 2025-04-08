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


def main(route):
    """
    This function detects if the synthesis follows a linear strategy.

    A linear synthesis strategy is characterized by:
    1. Multiple reaction steps (at least 3)
    2. No more than 2 non-trivial reactants per step
    3. A dominant pathway where most reactions have one main product and one main reactant
    """
    # Track reaction depths and branching
    reaction_depths = set()
    max_reactants_per_step = 0

    # Track the longest path
    path_lengths = []

    def dfs_traverse(node, depth=0):
        nonlocal reaction_depths, max_reactants_per_step

        if node["type"] == "reaction":
            # Store the depth
            reaction_depths.add(depth)

            # Extract reactants
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")

            # Count non-trivial reactants (excluding small molecules)
            non_trivial_reactants = 0
            for reactant in reactants_smiles:
                mol = Chem.MolFromSmiles(reactant)
                if mol and mol.GetNumHeavyAtoms() > 6:  # More reliable than SMILES length
                    non_trivial_reactants += 1

            max_reactants_per_step = max(max_reactants_per_step, non_trivial_reactants)

            # Process children and track path lengths
            child_paths = []
            for child in node.get("children", []):
                child_path = dfs_traverse(child, depth + 1)
                child_paths.append(child_path)

            # Return the longest path through this node
            if child_paths:
                return 1 + max(child_paths)
            return 1

        elif node["type"] == "mol":
            # For leaf nodes (starting materials)
            if not node.get("children", []):
                path_lengths.append(depth)
                return 0

            # For intermediate molecules, continue traversal
            child_paths = []
            for child in node.get("children", []):
                child_path = dfs_traverse(child, depth)
                child_paths.append(child_path)

            # Return the longest path through this node
            if child_paths:
                return max(child_paths)
            return 0

    # Start traversal and get the longest path
    longest_path = dfs_traverse(route)

    # Calculate the ratio of the longest path to total reactions
    path_ratio = longest_path / len(reaction_depths) if reaction_depths else 0

    # Linear synthesis typically has:
    # 1. Multiple reaction steps (at least 3)
    # 2. No more than 2 non-trivial reactants per step
    # 3. A high ratio of longest path to total reactions (indicating limited branching)
    result = (
        len(reaction_depths) >= 3 and max_reactants_per_step <= 2 and path_ratio >= 0.7
    )  # At least 70% of reactions are on the main path

    print(f"Linear synthesis strategy detected: {result}")
    print(f"Number of reaction steps: {len(reaction_depths)}")
    print(f"Max non-trivial reactants per step: {max_reactants_per_step}")
    print(f"Longest path ratio: {path_ratio:.2f}")

    return result
