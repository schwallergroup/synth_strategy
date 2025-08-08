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
    This function detects a linear synthesis strategy (as opposed to convergent).
    """
    # Track branching in the synthesis tree
    max_children = 0

    def dfs_traverse(node):
        nonlocal max_children

        if node["type"] == "reaction":
            # Count number of reactants (children)
            num_children = len(node.get("children", []))
            max_children = max(max_children, num_children)
            print(f"Reaction with {num_children} reactants found")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Linear synthesis typically has at most 2 reactants per step
    strategy_present = max_children <= 2

    print(f"Linear synthesis strategy detected: {strategy_present}")
    return strategy_present
