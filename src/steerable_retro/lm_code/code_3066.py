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
    Detects if the synthesis follows a linear (non-convergent) approach.

    A linear synthesis is characterized by:
    1. Each reaction step has at most one non-terminal reactant
    2. The overall structure forms a single main reaction pathway

    Returns:
    - True if the synthesis is convergent (non-linear)
    - False if the synthesis is linear
    """
    # Track if we've found any non-linear patterns
    is_non_linear = False

    def dfs_traverse(node, depth=0):
        nonlocal is_non_linear

        # For reaction nodes, check if there are multiple non-terminal reactants
        if node["type"] == "reaction":
            non_terminal_reactants = 0

            for child in node.get("children", []):
                # Skip terminal nodes (starting materials)
                if child.get("in_stock", False):
                    continue

                # Count non-terminal reactants (those that have their own synthesis path)
                if child["type"] == "mol" and not child.get("in_stock", False):
                    non_terminal_reactants += 1

                # If a reactant is itself a reaction, that's also non-linear
                elif child["type"] == "reaction":
                    non_terminal_reactants += 1

            # If we have more than one non-terminal reactant, it's a convergent (non-linear) synthesis
            if non_terminal_reactants > 1:
                is_non_linear = True
                print(
                    f"Non-linear pattern detected at depth {depth}: {non_terminal_reactants} non-terminal reactants"
                )

        # Continue traversing the tree
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Return True for convergent (non-linear) synthesis, False for linear synthesis
    is_convergent = is_non_linear

    if is_convergent:
        print("Detected convergent (non-linear) synthesis strategy")
    else:
        print("Detected linear synthesis strategy")

    return is_convergent
