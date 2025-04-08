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
    Detects late-stage Suzuki coupling connecting two heterocyclic fragments.
    Late stage is defined as occurring in the first half of the synthesis (lower depth).
    """
    max_depth = 0
    suzuki_depths = []

    # First pass to find max depth
    def find_max_depth(node, current_depth=0):
        nonlocal max_depth
        max_depth = max(max_depth, current_depth)

        for child in node.get("children", []):
            find_max_depth(child, current_depth + 1)

    # Second pass to find Suzuki couplings
    def find_suzuki_couplings(node, current_depth=0):
        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for boronic acid/ester pattern in reactants
            boronic_pattern = Chem.MolFromSmarts("[#6]-[#5](-[#8])-[#8]")

            for reactant in reactants:
                mol = Chem.MolFromSmiles(reactant)
                if mol and mol.HasSubstructMatch(boronic_pattern):
                    suzuki_depths.append(current_depth)
                    print(f"Found Suzuki coupling at depth {current_depth}")
                    break

        for child in node.get("children", []):
            find_suzuki_couplings(child, current_depth + 1)

    find_max_depth(route)
    find_suzuki_couplings(route)

    # Check if any Suzuki couplings occur in the first half of the synthesis
    for depth in suzuki_depths:
        if depth <= max_depth / 2:
            return True

    return False
