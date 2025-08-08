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


def main(route):
    """
    This function detects if the synthesis follows a linear strategy rather than convergent.
    Linear synthesis is characterized by having a single reactant in most reactions
    and a predominantly linear tree structure.
    """
    reaction_count = 0
    linear_reaction_count = 0

    # Track branching in the synthesis tree
    branching_nodes = 0
    total_nodes = 0

    def dfs_traverse(node, depth=0):
        nonlocal reaction_count, linear_reaction_count, branching_nodes, total_nodes

        total_nodes += 1

        # Count children for branching analysis
        child_count = len(node.get("children", []))

        if child_count > 1:
            branching_nodes += 1
            print(f"Found branching node at depth {depth} with {child_count} children")

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                reactants = reactants_part.split(".")

                reaction_count += 1

                # A reaction is considered linear if it has only one reactant
                if len(reactants) == 1:
                    linear_reaction_count += 1
                    print(f"Found linear reaction step: {rsmi}")
            except (IndexError, ValueError) as e:
                print(f"Error processing reaction SMILES: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Calculate linearity based on both reactant count and tree structure
    reactant_linearity = linear_reaction_count / reaction_count if reaction_count > 0 else 0
    structure_linearity = 1 - (branching_nodes / total_nodes) if total_nodes > 0 else 0

    # Combined linearity score (weighted average)
    linearity_score = (reactant_linearity + 2 * structure_linearity) / 3

    print(
        f"Total reactions: {reaction_count}, Linear reactions: {linear_reaction_count}, Reactant linearity: {reactant_linearity:.2f}"
    )
    print(
        f"Total nodes: {total_nodes}, Branching nodes: {branching_nodes}, Structure linearity: {structure_linearity:.2f}"
    )
    print(f"Overall linearity score: {linearity_score:.2f}")

    # Consider the synthesis linear if the combined score is above 0.5
    return linearity_score > 0.5
