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
    Detects if the synthesis follows a linear pattern without convergent steps.
    """
    # Track branching in the synthesis tree
    max_children_per_node = 0

    def dfs_traverse(node):
        nonlocal max_children_per_node

        if node["type"] == "reaction":
            # Count reactants
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            reactants = reactants_part.split(".")

            # Update max children count
            num_reactants = len([r for r in reactants if r])
            max_children_per_node = max(max_children_per_node, num_reactants)

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Linear synthesis has at most 1 non-reagent reactant per step
    if max_children_per_node <= 2:  # Allow for 1 main reactant + 1 reagent
        print("Linear synthesis pattern detected")
        return True

    print(f"Convergent synthesis detected with {max_children_per_node} maximum reactants")
    return False
