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
    This function detects if the synthesis uses a persistent Boc protecting group
    that remains on a nitrogen throughout multiple steps.
    """
    boc_protected_depths = []

    def dfs_traverse(node):
        if node["type"] == "reaction":
            depth = node.get("metadata", {}).get("depth", 0)
            rsmi = node.get("metadata", {}).get("rsmi", "")

            if not rsmi:
                return

            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for Boc group in product
            product_mol = Chem.MolFromSmiles(product)
            if product_mol:
                boc_pattern = Chem.MolFromSmarts(
                    "[N][C](=[O])[O][C]([CH3])([CH3])[CH3]"
                )
                if product_mol.HasSubstructMatch(boc_pattern):
                    boc_protected_depths.append(depth)

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # Check if Boc protection spans multiple consecutive steps
    if len(boc_protected_depths) >= 3:
        print(f"Found persistent Boc protection across depths: {boc_protected_depths}")
        return True
    return False
