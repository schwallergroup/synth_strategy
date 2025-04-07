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
    Detects if the synthesis uses a late-stage Suzuki coupling (depth 0 or 1)
    to connect two complex fragments.
    """
    suzuki_at_late_stage = False

    def dfs_traverse(node, depth=0):
        nonlocal suzuki_at_late_stage

        if node["type"] == "reaction" and depth <= 1:  # Late stage (depth 0 or 1)
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for boronic acid/ester pattern in reactants
            boronic_pattern = Chem.MolFromSmarts("[B;$(B-O)]")
            # Check for aryl halide pattern in reactants
            aryl_halide_pattern = Chem.MolFromSmarts("[c][Br,I,Cl]")

            has_boronic = False
            has_aryl_halide = False

            for reactant in reactants:
                mol = Chem.MolFromSmiles(reactant)
                if mol:
                    if mol.HasSubstructMatch(boronic_pattern):
                        has_boronic = True
                    if mol.HasSubstructMatch(aryl_halide_pattern):
                        has_aryl_halide = True

            # If both patterns are found, it's likely a Suzuki coupling
            if has_boronic and has_aryl_halide:
                print(f"Found late-stage Suzuki coupling at depth {depth}")
                suzuki_at_late_stage = True

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    return suzuki_at_late_stage
