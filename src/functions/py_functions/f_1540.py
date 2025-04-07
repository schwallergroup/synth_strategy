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
    Detects the formation of an O-alkylated hydroxylamine (oxime) structure
    in an early stage of the synthesis.
    """
    # Track oxime formation
    has_oxime_formation = False

    # SMARTS pattern for O-alkylated hydroxylamine (oxime)
    oxime_pattern = "[#6]-[#8]-[#7]=[#6]"

    def dfs_traverse(node, depth=0):
        nonlocal has_oxime_formation

        if node["type"] == "reaction" and depth >= 3:  # Early in synthesis
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                product = rsmi.split(">")[-1]

                product_mol = Chem.MolFromSmiles(product)
                if product_mol and product_mol.HasSubstructMatch(
                    Chem.MolFromSmarts(oxime_pattern)
                ):
                    has_oxime_formation = True
                    print(
                        f"Detected O-alkylated hydroxylamine formation at depth {depth}"
                    )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    if has_oxime_formation:
        print("Detected O-alkylated hydroxylamine formation strategy")

    return has_oxime_formation
