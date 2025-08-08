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
    Detects if the synthetic route follows a linear strategy (each intermediate molecule
    is the product of exactly one reaction).

    In a retrosynthetic tree:
    - The root is the target molecule
    - Each molecule node (except starting materials) should have exactly one reaction child
    - Reaction nodes can have multiple molecule children (reactants)
    """
    is_linear = True

    def dfs_traverse(node, depth=0):
        nonlocal is_linear

        if not is_linear:
            return  # Early return if we already know it's not linear

        if node["type"] == "mol" and not node.get("in_stock", False):
            # For non-starting molecules, check if they have exactly one reaction child
            reaction_children = [
                child for child in node.get("children", []) if child["type"] == "reaction"
            ]

            if len(reaction_children) != 1:
                is_linear = False
                print(
                    f"Non-linear synthesis detected (molecule with {len(reaction_children)} reactions)"
                )
                return

        # Continue traversing the tree
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    if is_linear:
        print("Linear synthesis strategy detected")

    return is_linear
