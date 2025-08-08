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
    This function detects if the synthesis follows a linear pattern
    rather than a convergent approach.
    """
    # Track number of reactants per step
    reactant_counts = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            # Extract reactants
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")

            # Count non-empty reactants
            count = sum(1 for r in reactants if r.strip())
            reactant_counts.append(count)
            print(f"Reaction at depth {depth} has {count} reactants")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # If most reactions have 2 reactants, it's likely a linear synthesis
    # (one main substrate + one reagent in each step)
    if (
        reactant_counts
        and sum(count == 2 for count in reactant_counts) / len(reactant_counts) >= 0.7
    ):
        print("Detected linear synthesis pattern")
        return True
    return False
