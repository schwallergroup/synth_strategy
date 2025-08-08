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
    This function detects a linear build-up strategy where complexity is added
    sequentially rather than through convergent synthesis.
    """
    reaction_depths = []
    max_reactants_per_step = 0

    def dfs_traverse(node, depth=0):
        nonlocal reaction_depths, max_reactants_per_step

        if node["type"] == "reaction":
            reaction_depths.append(depth)

            # Extract reactants
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")

            # Count number of reactants
            num_reactants = len(reactants_smiles)
            max_reactants_per_step = max(max_reactants_per_step, num_reactants)

        # Continue traversing children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if we have a linear build-up (consecutive depths, limited reactants per step)
    is_linear = len(reaction_depths) >= 3 and max_reactants_per_step <= 2

    if is_linear:
        print("Detected linear build-up strategy")

    return is_linear
