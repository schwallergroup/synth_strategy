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
    This function detects if the synthetic route follows a linear synthesis strategy
    (each step builds on a single previous product).
    """
    is_linear = True
    reaction_count = 0

    def dfs_traverse(node):
        nonlocal is_linear, reaction_count

        if node["type"] == "reaction":
            reaction_count += 1
            # Extract reactants
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")

            # If there are more than 2 reactants, it's likely not a linear synthesis
            # (allowing for up to 2 to account for reagents)
            if len(reactants_smiles) > 2:
                is_linear = False
                print(
                    f"Found non-linear step at depth {node.get('depth', 'unknown')} with {len(reactants_smiles)} reactants"
                )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Only consider routes with multiple reactions
    result = is_linear and reaction_count > 1
    print(f"Linear synthesis strategy detected: {result} (reaction count: {reaction_count})")
    return result
