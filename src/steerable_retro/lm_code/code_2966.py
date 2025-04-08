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
    Detects if the synthesis uses a late-stage Sonogashira coupling (C-C≡C bond formation)
    between two complex fragments.
    """
    sonogashira_detected = False
    late_stage = False

    def dfs_traverse(node, depth=0):
        nonlocal sonogashira_detected, late_stage

        if node["type"] == "reaction" and depth <= 1:  # Late stage (depth 0 or 1)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if this is potentially a Sonogashira coupling
                # Look for C≡C bond in product that connects two fragments
                if "#" in product:
                    # Check if we're combining two fragments
                    if len(reactants) >= 2:
                        # Check if one fragment has a halide (I, Br) and another has an alkyne
                        has_halide = False
                        has_alkyne = False

                        for reactant in reactants:
                            if "I" in reactant or "Br" in reactant:
                                has_halide = True
                            if "#" in reactant:
                                has_alkyne = True

                        if has_halide and has_alkyne:
                            sonogashira_detected = True
                            late_stage = True
                            print(f"Late-stage Sonogashira coupling detected at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return sonogashira_detected and late_stage
