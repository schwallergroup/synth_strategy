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
    This function detects a late-stage Suzuki cross-coupling strategy where
    boronic ester and triflate/halide coupling partners are prepared and then joined.
    """
    # Track if we found the key components
    found_borylation = False
    found_triflation_or_halide = False
    found_coupling = False
    reaction_depths = {}
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal found_borylation, found_triflation_or_halide, found_coupling, max_depth

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            reaction_depths[depth] = reaction_depths.get(depth, 0) + 1

            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for borylation (aryl-Br → aryl-B(OR)2)
                if any("Br" in r for r in reactants) and "B" in product and "O" in product:
                    found_borylation = True
                    print(f"Found borylation at depth {depth}")

                # Check for triflation or halide formation
                if any("S(=O)(=O)" in r and "F" in r for r in reactants) or any(
                    "Br" in product for r in reactants
                ):
                    found_triflation_or_halide = True
                    print(f"Found triflation/halide at depth {depth}")

                # Check for coupling (two fragments joining)
                if (
                    len(reactants) >= 2
                    and "B" in rsmi
                    and ("O" in rsmi or "Br" in rsmi or "I" in rsmi)
                ):
                    # This is a simplification - real implementation would check for actual coupling
                    found_coupling = True
                    print(f"Found potential coupling reaction at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if coupling occurs in the second half of synthesis
    late_stage = False
    if found_coupling:
        for depth, count in reaction_depths.items():
            if depth <= max_depth / 2:  # First half of synthesis (lower depth values)
                late_stage = True
                break

    result = found_borylation and found_triflation_or_halide and found_coupling and late_stage
    print(f"Late-stage Suzuki coupling strategy detected: {result}")
    return result
