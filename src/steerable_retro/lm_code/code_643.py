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
    This function detects if the synthesis follows a linear strategy rather than convergent.
    Linear synthesis typically has 1-2 reactants per step.
    """
    is_linear = True
    step_count = 0

    def dfs_traverse(node, depth=0):
        nonlocal is_linear, step_count

        if node["type"] == "reaction":
            step_count += 1
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")

            # Count non-empty reactants
            reactant_count = sum(1 for r in reactants if r.strip())

            # If more than 2 significant reactants, it's likely a convergent step
            if reactant_count > 2:
                print(f"Convergent step detected at depth {depth} with {reactant_count} reactants")
                is_linear = False

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Need at least 3 steps to meaningfully classify as linear
    if step_count < 3:
        print(f"Not enough steps ({step_count}) to classify synthesis strategy")
        return False

    if is_linear:
        print("Linear synthesis strategy detected")

    return is_linear
