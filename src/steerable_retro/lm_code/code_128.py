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
    """Helper function to count the total number of reaction nodes"""
    reaction_count = 0

    def dfs_count(node):
        nonlocal reaction_count
        if node["type"] == "reaction":
            reaction_count += 1
        for child in node.get("children", []):
            dfs_count(child)

    dfs_count(route)
    return reaction_count
