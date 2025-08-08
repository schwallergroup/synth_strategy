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
    This function detects a linear synthesis strategy where each reaction
    typically has one main building block added to a growing molecule.
    """
    linear_steps = 0
    total_steps = 0

    def dfs_traverse(node, depth=0):
        nonlocal linear_steps, total_steps

        if node["type"] == "reaction":
            total_steps += 1

            # Check number of reactants
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # Linear steps typically have 1-2 reactants
                if len(reactants) <= 2:
                    linear_steps += 1
                    print(f"Linear step detected at depth {depth} with {len(reactants)} reactants")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if most steps are linear (at least 75%)
    if total_steps > 0 and linear_steps / total_steps >= 0.75:
        print(f"Linear synthesis strategy detected: {linear_steps}/{total_steps} steps are linear")
        return True
    return False
