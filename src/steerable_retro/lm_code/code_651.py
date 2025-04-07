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
    This function detects if a morpholine group is introduced in the second half of the synthesis.
    """
    morpholine_introduction_detected = False
    max_depth = 0

    # First pass to determine the maximum depth
    def get_max_depth(node, depth=0):
        nonlocal max_depth
        max_depth = max(max_depth, depth)
        for child in node.get("children", []):
            get_max_depth(child, depth + 1)

    get_max_depth(route)

    # Second pass to detect morpholine introduction
    def dfs_traverse(node, depth=0):
        nonlocal morpholine_introduction_detected

        if node["type"] == "reaction" and depth < max_depth / 2:  # Second half of synthesis
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for morpholine in reactants
                morpholine_pattern = Chem.MolFromSmarts("[N]1CCO[C][C]1")
                morpholine_in_reactants = False
                for reactant in reactants:
                    if (
                        reactant
                        and Chem.MolFromSmiles(reactant)
                        and Chem.MolFromSmiles(reactant).HasSubstructMatch(morpholine_pattern)
                    ):
                        morpholine_in_reactants = True
                        break

                # Check for morpholine in product but not in all reactants
                morpholine_in_product = False
                if (
                    product
                    and Chem.MolFromSmiles(product)
                    and Chem.MolFromSmiles(product).HasSubstructMatch(morpholine_pattern)
                ):
                    morpholine_in_product = True

                if morpholine_in_reactants and morpholine_in_product:
                    print("Morpholine introduction detected at depth", depth)
                    morpholine_introduction_detected = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Call dfs_traverse on the root node
    dfs_traverse(route)

    return morpholine_introduction_detected
