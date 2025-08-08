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
    This function detects if the synthesis follows a linear strategy with no branching.

    In a linear synthesis:
    - Each molecule node should lead to at most one reaction node
    - Reaction nodes can have multiple children (reactants) without affecting linearity
    """
    is_linear = True

    def dfs_traverse(node):
        nonlocal is_linear

        # Only check branching for molecule nodes
        if node.get("type") == "mol" and not node.get("in_stock", False):
            reaction_children = [
                child for child in node.get("children", []) if child.get("type") == "reaction"
            ]

            # If a molecule leads to more than one reaction, it's not linear
            if len(reaction_children) > 1:
                print(f"Found branching at molecule: {node.get('smiles', 'unknown')}")
                is_linear = False

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    print(f"Linear synthesis strategy: {is_linear}")
    return is_linear
