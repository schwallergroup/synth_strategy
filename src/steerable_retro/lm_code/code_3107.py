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
    Detects if the synthesis follows a linear pattern without convergent branches
    """
    # Track the maximum branching factor
    max_branching = 0

    def dfs_traverse(node, depth=0):
        nonlocal max_branching

        if node["type"] == "reaction":
            # Count the number of reactants
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            branching = len(reactants)
            max_branching = max(max_branching, branching)

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # If max_branching is consistently 2 or less, it's likely a linear synthesis
    # (one main reactant plus one reagent/functional group)
    is_linear = max_branching <= 2
    print(f"Maximum branching factor: {max_branching}, Linear synthesis: {is_linear}")

    return is_linear
