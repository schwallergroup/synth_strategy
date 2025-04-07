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
    Detects multiple Suzuki coupling reactions in the synthetic route.
    Looks for C-C bond formations with boronic acid reactants.
    """
    suzuki_coupling_depths = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")

            # Check for boronic acid pattern in reactants
            if any("OB(O)" in r for r in reactants):
                suzuki_coupling_depths.append(depth)
                print(f"Found Suzuki coupling at depth {depth}: {rsmi}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    if len(suzuki_coupling_depths) >= 2:
        print(f"Multiple Suzuki couplings detected at depths: {suzuki_coupling_depths}")
        return True
    return False
