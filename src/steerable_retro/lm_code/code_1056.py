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
    Detects if the synthetic route follows a linear fragment assembly pattern
    with at least 3 distinct fragments
    """
    fragment_count = 0
    is_linear = True

    def dfs_traverse(node, depth=0):
        nonlocal fragment_count, is_linear

        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # Count number of reactants as a proxy for fragments
                if len(reactants) > 1:
                    fragment_count += 1

                # If any reaction has more than 2 reactants, it might not be linear
                if len(reactants) > 2:
                    is_linear = False

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # A linear assembly should have at least 3 fragment combinations
    # and maintain a linear pattern (mostly 2 reactants per step)
    result = fragment_count >= 3 and is_linear
    if result:
        print(f"Found linear fragment assembly with {fragment_count} fragment combinations")

    return result
