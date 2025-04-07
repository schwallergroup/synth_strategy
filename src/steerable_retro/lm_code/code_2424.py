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
    Detects a synthetic strategy involving late-stage Suzuki coupling (depth 0-1)
    with boronic acid/ester and aryl halide.
    """
    suzuki_coupling_found = False

    def dfs_traverse(node, depth=0):
        nonlocal suzuki_coupling_found

        if node["type"] == "reaction" and depth <= 1:  # Late stage (depth 0-1)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # Check for boronic acid/ester pattern in reactants
                boronic_pattern = re.compile(r"B.*O", re.IGNORECASE)
                has_boronic = any(boronic_pattern.search(r) for r in reactants)

                # Check for aryl halide pattern in reactants
                aryl_halide_pattern = re.compile(r"c.*[BrClI]", re.IGNORECASE)
                has_aryl_halide = any(aryl_halide_pattern.search(r) for r in reactants)

                if has_boronic and has_aryl_halide:
                    print("Found Suzuki coupling at depth", depth)
                    suzuki_coupling_found = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return suzuki_coupling_found
