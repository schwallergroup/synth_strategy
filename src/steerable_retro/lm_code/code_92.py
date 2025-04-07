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
    This function detects if the synthesis route employs metal-mediated coupling reactions,
    specifically looking for organozinc reagents.
    """
    has_metal_coupling = False

    def dfs_traverse(node, depth=0):
        nonlocal has_metal_coupling

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")

            # Check for zinc reagents in reactants
            for reactant in reactants:
                if "[Zn+]" in reactant:
                    has_metal_coupling = True
                    print(f"Found metal-mediated coupling with zinc at depth {depth}")
                    break

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    print(f"Has metal-mediated coupling: {has_metal_coupling}")
    return has_metal_coupling
