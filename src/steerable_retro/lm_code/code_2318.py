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
    Detects if the synthesis route uses a late-stage Suzuki coupling (depth 0 or 1)
    between a boronic acid/ester and an aryl halide.
    """
    suzuki_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal suzuki_detected

        if node["type"] == "reaction" and depth <= 1:  # Only check late-stage reactions
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for boronic acid/ester pattern in reactants
                boronic_pattern = Chem.MolFromSmarts("[#6]-[#5](-[#8])-[#8]")
                # Check for aryl halide pattern in reactants
                aryl_halide_pattern = Chem.MolFromSmarts("[#6]-[#17,#35,#53]")

                has_boronic = False
                has_aryl_halide = False

                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol:
                        if mol.HasSubstructMatch(boronic_pattern):
                            has_boronic = True
                        if mol.HasSubstructMatch(aryl_halide_pattern):
                            has_aryl_halide = True

                if has_boronic and has_aryl_halide:
                    print("Detected Suzuki coupling at depth", depth)
                    suzuki_detected = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    return suzuki_detected
