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
    Detects if the synthesis route follows a predominantly linear strategy.
    A linear synthesis has most reaction nodes with only one non-commercial reactant.
    """
    reaction_nodes = []

    def collect_reaction_nodes(node, depth=0):
        if node["type"] == "reaction":
            reaction_nodes.append((node, depth))
        for child in node.get("children", []):
            collect_reaction_nodes(child, depth + 1)

    collect_reaction_nodes(route)

    if not reaction_nodes:
        return False

    linear_count = 0
    for node, _ in reaction_nodes:
        non_commercial_reactants = 0
        for child in node.get("children", []):
            if child["type"] == "mol" and not child.get("in_stock", False):
                non_commercial_reactants += 1

        if non_commercial_reactants <= 1:
            linear_count += 1

    # If at least 70% of reactions have at most one non-commercial reactant, consider it linear
    linearity_ratio = linear_count / len(reaction_nodes)
    print(
        f"Linearity ratio: {linearity_ratio:.2f} ({linear_count}/{len(reaction_nodes)} reactions)"
    )

    return linearity_ratio >= 0.7
