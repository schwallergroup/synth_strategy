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
    This function detects if the synthesis follows a linear (non-convergent) strategy.
    """
    # Track the maximum branching factor in the reaction tree
    max_branching = 0

    def dfs_traverse(node):
        nonlocal max_branching

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # Count number of reactants (branching factor)
                branching = len(reactants)
                max_branching = max(max_branching, branching)

                if branching > 1:
                    print(f"Found reaction with {branching} reactants")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from root
    dfs_traverse(route)

    # If max branching is 1, it's a strictly linear synthesis
    # If max branching is 2, it's mostly linear with some convergent steps
    # If max branching > 2, it's a more convergent synthesis
    is_linear = max_branching <= 2
    print(f"Maximum branching factor: {max_branching}, Linear synthesis: {is_linear}")
    return is_linear
