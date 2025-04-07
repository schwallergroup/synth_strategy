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
    Detects late-stage functional group transformation (depth 0 or 1)
    """
    late_stage_transformation = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_transformation

        if node["type"] == "reaction" and depth <= 1:
            # This is a late-stage transformation (depth 0 or 1)
            print(f"Late-stage transformation detected at depth {depth}")
            late_stage_transformation = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return late_stage_transformation
