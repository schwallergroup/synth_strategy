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
    This function detects if the synthesis follows a linear strategy (each step builds on the previous
    without convergent steps).
    """
    is_linear = True

    def count_children(node):
        if node["type"] == "reaction":
            children = node.get("children", [])
            # If a reaction has more than 2 children (more than 1 non-trivial reactant),
            # it's likely a convergent step
            non_trivial_children = 0
            for child in children:
                if child["type"] == "mol" and not child.get("in_stock", False):
                    non_trivial_children += 1
                elif child["type"] == "reaction":
                    non_trivial_children += 1

            if non_trivial_children > 1:
                return False

            # Recursively check children
            for child in children:
                if not count_children(child):
                    return False

        return True

    is_linear = count_children(route)
    if is_linear:
        print("Linear synthesis strategy detected")
    else:
        print("Convergent synthesis strategy detected")

    return is_linear
