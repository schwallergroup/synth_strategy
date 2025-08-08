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
    Detects if the synthesis follows a linear pattern (no convergent steps).
    """
    # Track the maximum branching factor in the synthesis tree
    max_branching = 0

    def dfs_traverse(node):
        nonlocal max_branching

        if node["type"] == "reaction":
            # Count number of reactants (children that are molecules)
            reactant_count = sum(1 for child in node.get("children", []) if child["type"] == "mol")
            max_branching = max(max_branching, reactant_count)

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # If max_branching <= 2, it's likely a linear synthesis
    # (allowing for 2 to account for reagents that aren't part of the main skeleton)
    is_linear = max_branching <= 2
    print(f"Linear synthesis pattern: {is_linear} (max branching: {max_branching})")
    return is_linear
