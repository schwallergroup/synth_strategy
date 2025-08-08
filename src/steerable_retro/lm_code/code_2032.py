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
    This function detects if the synthesis follows a linear build-up approach
    rather than a convergent approach.
    """
    # Track the maximum number of reactants in any step
    max_reactants = 0
    # Track total number of steps
    total_steps = 0

    def dfs_traverse(node):
        nonlocal max_reactants, total_steps

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # Count number of reactants
                num_reactants = len(reactants)
                max_reactants = max(max_reactants, num_reactants)
                total_steps += 1

                print(f"Step with {num_reactants} reactants: {rsmi}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # A linear synthesis typically has fewer reactants per step (usually 1-2)
    # We'll define linear as having at most 2 reactants in any step and at least 3 steps
    is_linear = max_reactants <= 2 and total_steps >= 3

    print(
        f"Linear synthesis strategy detected: {is_linear} (max_reactants={max_reactants}, total_steps={total_steps})"
    )
    return is_linear
