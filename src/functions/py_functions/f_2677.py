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
    Detects if the synthesis follows a linear strategy where each reaction
    has only one product that feeds into the next reaction.

    In retrosynthetic analysis, we check if each reaction has only one non-stock
    reactant that continues the synthesis chain, and if each non-stock molecule
    feeds into only one reaction.
    """
    is_linear = True
    reaction_count = 0

    def dfs_traverse(node, depth=0):
        nonlocal is_linear, reaction_count

        if node["type"] == "reaction":
            # In retrosynthesis, children of a reaction are reactants
            non_stock_reactants = [
                child
                for child in node.get("children", [])
                if child["type"] == "mol" and not child.get("in_stock", False)
            ]

            # If more than one non-stock reactant, the synthesis is not linear
            if len(non_stock_reactants) > 1:
                is_linear = False
                print(
                    f"Non-linear step detected at depth {depth}: reaction has multiple non-stock reactants"
                )

            # Count this reaction if it's part of the main synthesis path
            if non_stock_reactants:
                reaction_count += 1

        elif node["type"] == "mol" and not node.get("in_stock", False):
            # If this molecule feeds into multiple reactions, the synthesis branches
            reaction_children = [
                child
                for child in node.get("children", [])
                if child["type"] == "reaction"
            ]

            if len(reaction_children) > 1:
                is_linear = False
                print(
                    f"Non-linear step detected at depth {depth}: molecule branches into {len(reaction_children)} reactions"
                )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # A linear synthesis should have at least 2 reactions
    if reaction_count < 2:
        print(f"Not enough reactions: found {reaction_count}, need at least 2")
        return False

    return is_linear
