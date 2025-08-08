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
    This function detects a linear synthesis strategy where each step builds on a single precursor
    rather than combining multiple complex fragments.
    """
    # Initialize tracking variables
    is_linear = True
    max_reactants_per_step = 0

    def dfs_traverse(node, depth=0):
        nonlocal is_linear, max_reactants_per_step

        if node["type"] == "reaction":
            # Extract reactants
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")

            # Count non-empty reactants
            num_reactants = sum(1 for r in reactants_smiles if r)
            max_reactants_per_step = max(max_reactants_per_step, num_reactants)

            # If any step has more than 2 reactants, it's likely not a linear synthesis
            if num_reactants > 2:
                is_linear = False
                print(f"Non-linear step detected at depth {depth} with {num_reactants} reactants")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # A truly linear synthesis typically has at most 2 reactants per step
    # (one main building block and one reagent/coupling partner)
    is_linear = is_linear and max_reactants_per_step <= 2

    print(f"Maximum reactants per step: {max_reactants_per_step}")
    print(f"Linear synthesis strategy present: {is_linear}")

    return is_linear
