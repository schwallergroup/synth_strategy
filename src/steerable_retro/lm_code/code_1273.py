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
    Linear synthesis is characterized by each reaction having only one non-commercial reactant.
    """
    is_linear = True

    def dfs_traverse(node):
        nonlocal is_linear

        if node["type"] == "reaction":
            # Count non-commercial (non-leaf) reactants
            non_commercial_reactants = 0

            for child in node.get("children", []):
                if child["type"] == "mol":
                    if not child.get("in_stock", False):
                        non_commercial_reactants += 1

            # If more than one non-commercial reactant, it's not strictly linear
            if non_commercial_reactants > 1:
                is_linear = False
                print(
                    f"Non-linear (convergent) step detected with {non_commercial_reactants} non-commercial reactants"
                )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Call dfs_traverse on the root node
    dfs_traverse(route)

    if is_linear:
        print("Synthesis follows a linear strategy")

    return is_linear
