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
    This function detects if the synthesis follows a predominantly linear strategy
    (as opposed to convergent) by checking if most reactions have only one non-reagent reactant.
    """
    reaction_count = 0
    linear_reaction_count = 0

    def dfs_traverse(node):
        nonlocal reaction_count, linear_reaction_count

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")

            reaction_count += 1

            # Count significant fragments (excluding small reagents)
            significant_fragments = 0
            for reactant in reactants:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if (
                        mol and mol.GetNumHeavyAtoms() > 6
                    ):  # Arbitrary threshold for "significant" fragment
                        significant_fragments += 1
                except:
                    pass

            if significant_fragments <= 1:
                linear_reaction_count += 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    # Consider it a linear strategy if most reactions are linear
    if reaction_count > 0 and linear_reaction_count / reaction_count > 0.7:
        print(
            f"Linear synthesis strategy detected ({linear_reaction_count}/{reaction_count} reactions)"
        )
        return True
    return False
