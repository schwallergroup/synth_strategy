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
    This function detects a linear build-up strategy where complexity
    is added sequentially without convergent steps.
    """
    reaction_count = 0
    max_reactants_per_step = 0

    def dfs_traverse(node, depth=0):
        nonlocal reaction_count, max_reactants_per_step

        if node["type"] == "reaction":
            reaction_count += 1
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            max_reactants_per_step = max(max_reactants_per_step, len(reactants))

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # If we have multiple reactions but never more than 2 reactants per step,
    # it's likely a linear build-up strategy
    if reaction_count >= 3 and max_reactants_per_step <= 2:
        print(f"Linear build-up strategy detected with {reaction_count} reactions")
        print(f"Maximum reactants per step: {max_reactants_per_step}")
        return True

    return False
