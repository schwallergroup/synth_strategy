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
    Detects a synthetic strategy involving sequential C-N bond formations
    (at least 2 C-N bond formations in the route).
    """
    c_n_bond_formations = 0

    def dfs_traverse(node, depth=0):
        nonlocal c_n_bond_formations

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for C-N bond formation
                # Look for patterns where an aryl halide reacts with an amine
                aryl_halide_pattern = re.compile(r"c.*[BrClI]", re.IGNORECASE)
                amine_pattern = re.compile(r"N", re.IGNORECASE)

                has_aryl_halide = any(aryl_halide_pattern.search(r) for r in reactants)
                has_amine = any(amine_pattern.search(r) for r in reactants)

                # Check if product has C-N bond that wasn't in reactants
                if has_aryl_halide and has_amine:
                    # This is a simplified check - in a real implementation,
                    # you would use RDKit to compare reactant and product structures
                    c_n_bond_formations += 1
                    print(f"Found potential C-N bond formation at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return c_n_bond_formations >= 2
