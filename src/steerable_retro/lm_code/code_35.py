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
    In a linear synthesis, each reaction has only one product that serves as a reactant in the next step.
    """
    # Track reaction depths and number of reactants at each depth
    reaction_depths = {}

    def dfs_traverse(node, depth=0):
        if node.get("type") == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")

            # Store the number of reactants at this depth
            reaction_depths[depth] = len(reactants)

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if most reactions have only one product feeding into the next step
    # For a linear synthesis, most depths should have only one reactant
    linear_steps = sum(1 for count in reaction_depths.values() if count <= 2)
    convergent_steps = sum(1 for count in reaction_depths.values() if count > 2)

    is_linear = linear_steps > convergent_steps

    if is_linear:
        print("Detected linear synthesis strategy")
    else:
        print("Detected convergent synthesis strategy")

    return is_linear
