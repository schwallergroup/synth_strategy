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
    Detects if the synthesis preserves a dimethoxyphenyl group throughout the route.
    """
    dimethoxyphenyl_pattern = Chem.MolFromSmarts("[c]1[c]([O][C])[c][c][c][c]1[O][C]")
    has_dimethoxyphenyl_at_all_depths = True
    depths_checked = set()

    def dfs_traverse(node):
        nonlocal has_dimethoxyphenyl_at_all_depths, depths_checked

        if node["type"] == "mol" and "smiles" in node:
            mol = Chem.MolFromSmiles(node["smiles"])
            depth = None

            # Try to get depth from parent reaction if this is a product
            if node.get("children"):
                for child in node["children"]:
                    if (
                        child["type"] == "reaction"
                        and "metadata" in child
                        and "depth" in child["metadata"]
                    ):
                        depth = child["metadata"]["depth"]
                        break

            if mol and depth is not None:
                depths_checked.add(depth)
                if not mol.HasSubstructMatch(dimethoxyphenyl_pattern):
                    has_dimethoxyphenyl_at_all_depths = False
                    print(f"Dimethoxyphenyl group not found at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    # Only return True if we've checked multiple depths and found the group at all of them
    return has_dimethoxyphenyl_at_all_depths and len(depths_checked) >= 3
