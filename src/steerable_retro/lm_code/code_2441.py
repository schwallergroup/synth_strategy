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
    This function detects a convergent synthesis approach with multiple fragments.
    """
    late_stage_fragment_combination = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_fragment_combination

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")

            # Check if this is a late-stage reaction (depth < 2) with multiple fragments
            if depth < 2 and len(reactants) >= 2:
                # Check if reactants are complex (more than 15 atoms)
                complex_reactants = sum(1 for r in reactants if len(r) > 15)
                if complex_reactants >= 2:
                    print(f"Found late-stage fragment combination at depth {depth}: {rsmi}")
                    late_stage_fragment_combination = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    return late_stage_fragment_combination
