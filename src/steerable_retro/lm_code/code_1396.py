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
    Detects if the synthesis follows a linear strategy without convergent steps
    """
    # Track the maximum branching factor at each depth
    branching_factors = {}

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            # Count number of reactants
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            num_reactants = len(reactants)

            # Update maximum branching factor at this depth
            if depth not in branching_factors or num_reactants > branching_factors[depth]:
                branching_factors[depth] = num_reactants

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    print(f"Branching factors at each depth: {branching_factors}")

    # A linear synthesis should have a maximum branching factor of 2
    # (allowing for reagents that aren't part of the main synthetic path)
    is_linear = all(bf <= 2 for bf in branching_factors.values())

    return is_linear
