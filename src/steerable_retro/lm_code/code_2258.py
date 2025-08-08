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
    """Check if the synthesis follows a predominantly linear strategy"""
    # Count branching factors at each node
    branching_factors = []

    def dfs(node):
        if node["type"] == "reaction":
            # Count number of reactants (children)
            num_children = len(node.get("children", []))
            if num_children > 0:
                branching_factors.append(num_children)

        # Recursively check children
        for child in node.get("children", []):
            dfs(child)

    dfs(route)

    # Calculate average branching factor
    if not branching_factors:
        return False

    avg_branching = sum(branching_factors) / len(branching_factors)
    is_linear = avg_branching <= 1.7  # Allow more branching

    print(f"Average branching factor: {avg_branching}, Linear synthesis: {is_linear}")
    return is_linear
