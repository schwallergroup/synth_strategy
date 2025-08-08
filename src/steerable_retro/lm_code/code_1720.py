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
    Detects if the route follows a linear synthesis strategy where each step
    builds directly on the previous product.
    """
    branching_factors = []

    def dfs_traverse(node):
        if node["type"] == "reaction":
            # Count number of reactants in this reaction
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                branching_factors.append(len(reactants))
                print(f"Reaction has {len(reactants)} reactants")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # If most reactions have only 1-2 reactants, it's likely a linear synthesis
    if branching_factors and sum(1 for bf in branching_factors if bf <= 2) >= 0.7 * len(
        branching_factors
    ):
        print("Linear synthesis strategy detected")
        return True
    return False
