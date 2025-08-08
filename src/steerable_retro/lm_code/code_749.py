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
    This function detects if the synthesis follows a linear rather than convergent approach.
    """
    # Track the maximum branching factor in the synthesis tree
    max_branching = 0

    def dfs_traverse(node):
        nonlocal max_branching

        if node["type"] == "reaction":
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if rsmi:
                reactants = rsmi.split(">")[0].split(".")
                num_reactants = len(reactants)
                max_branching = max(max_branching, num_reactants)

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # If max branching is 2 or less, consider it a linear synthesis
    is_linear = max_branching <= 2
    if is_linear:
        print("Linear synthesis strategy detected (max branching ≤ 2)")
    else:
        print(f"Convergent synthesis detected with max branching of {max_branching}")

    return is_linear
