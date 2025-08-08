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
    This function detects a linear fragment assembly strategy where fragments
    are sequentially added rather than converging multiple complex fragments.
    """
    reaction_depths = []

    def dfs_traverse(node, depth=0):
        nonlocal reaction_depths

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")

            # Store number of reactants at each depth
            reaction_depths.append((depth, len(reactants)))

        # Recursively traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Sort by depth
    reaction_depths.sort(key=lambda x: x[0])

    # Check if most reactions have 2 reactants (binary reactions)
    # and if they occur in sequence (linear assembly)
    binary_reactions = sum(1 for _, num_reactants in reaction_depths if num_reactants == 2)

    # Linear strategy typically has most reactions as binary and in sequence
    is_linear = binary_reactions >= len(reaction_depths) * 0.7

    if is_linear:
        print("Found linear fragment assembly strategy")

    return is_linear
