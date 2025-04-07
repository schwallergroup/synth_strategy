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
    Detects if the synthetic route contains multiple Suzuki coupling reactions.
    Suzuki coupling is identified by C-C bond formation between aryl groups where
    reactants contain boronic acid and aryl halide.
    """
    suzuki_count = 0

    def dfs_traverse(node):
        nonlocal suzuki_count

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for boronic acid in reactants
                has_boronic_acid = any(
                    "B(O)(O)" in r or "OB(O)" in r for r in reactants
                )
                # Check for aryl halide in reactants
                has_aryl_halide = any("Br" in r for r in reactants)

                if has_boronic_acid and has_aryl_halide:
                    suzuki_count += 1
                    print(f"Found Suzuki coupling: {rsmi}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(f"Total Suzuki couplings found: {suzuki_count}")
    return suzuki_count >= 2
