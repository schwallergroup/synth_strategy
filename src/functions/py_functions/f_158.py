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
    This function detects if the synthetic route involves late-stage aromatic bromination
    (in the second half of the synthesis).
    """
    bromination_detected = False
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal bromination_detected, max_depth

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this is a bromination reaction (product has bromine on aromatic carbon that wasn't there in reactants)
            product_mol = Chem.MolFromSmiles(product)
            if product_mol:
                # Look for aromatic carbon with bromine
                aromatic_br_pattern = Chem.MolFromSmarts("c[Br]")
                if product_mol.HasSubstructMatch(aromatic_br_pattern):
                    # Check if this is in the second half of synthesis (depth < max_depth/2)
                    if depth <= max_depth / 2:
                        print(
                            f"Late-stage aromatic bromination detected at depth {depth}"
                        )
                        bromination_detected = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # First pass to determine max_depth
    dfs_traverse(route)
    # Reset and do actual detection
    bromination_detected = False
    dfs_traverse(route)

    return bromination_detected
