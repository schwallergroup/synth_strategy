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
    This function detects if the synthetic route involves late-stage protection (phthalimide).
    """
    protection_depth = None
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal protection_depth, max_depth

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                product = rsmi.split(">")[-1]

                # Check for phthalimide in product
                phthalimide_pattern = Chem.MolFromSmarts("O=C1c2ccccc2C(=O)N1[C]")

                try:
                    prod_mol = Chem.MolFromSmiles(product)
                    if prod_mol and prod_mol.HasSubstructMatch(phthalimide_pattern):
                        protection_depth = depth
                        print(f"Phthalimide protection detected at depth {depth}")
                except:
                    pass

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # If protection occurs in the first half of the synthesis (lower depth numbers)
    if protection_depth is not None and protection_depth <= max_depth / 2:
        print(
            f"Late-stage protection detected (depth {protection_depth} of max {max_depth})"
        )
        return True
    return False
