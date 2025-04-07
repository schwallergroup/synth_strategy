#!/bin/python

"""LM-defined function for strategy description."""

import copy
import re
from collections import deque

import rdkit
import rdkit.Chem as Chem
from rdkit import Chem
from rdkit.Chem import (
    AllChem,
    Descriptors,
    Lipinski,
    rdChemReactions,
    rdFMCS,
    rdMolDescriptors,
    rdmolops,
)
from rdkit.Chem.Scaffolds import MurckoScaffold


def main(route):
    """
    Detects if the synthesis follows a linear strategy rather than convergent.
    """
    # Track the maximum branching factor in the synthesis tree
    max_branching = 0

    def dfs_traverse(node):
        nonlocal max_branching

        if node["type"] == "reaction":
            # Count number of reactants in this reaction
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0]
                reactants = reactants_smiles.split(".")
                branching = len(reactants)
                max_branching = max(max_branching, branching)

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    # If max branching is <= 2, consider it a linear synthesis
    is_linear = max_branching <= 2
    if is_linear:
        print("Linear synthesis strategy detected")
    else:
        print(f"Convergent synthesis detected with branching factor {max_branching}")

    return is_linear
