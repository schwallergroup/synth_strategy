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
    This function detects if the synthetic route employs a linear synthesis strategy
    rather than a convergent one. In linear synthesis, each step has one main reactant
    derived from the previous step.
    """
    # Count branches in the synthesis tree
    branch_counts = []

    def count_branches(node):
        if node["type"] == "reaction":
            children = node.get("children", [])
            branch_counts.append(len(children))

        for child in node.get("children", []):
            count_branches(child)

    count_branches(route)

    # If most reactions have only one or two reactants, it's likely a linear synthesis
    if branch_counts and sum(1 for c in branch_counts if c <= 2) / len(branch_counts) > 0.7:
        print("Detected linear synthesis strategy")
        return True
    return False
