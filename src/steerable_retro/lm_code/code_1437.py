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
    This function detects if the synthesis involves a Boc-protected amine
    that is maintained throughout multiple steps.
    """
    # Track Boc group presence at different depths
    boc_at_depths = set()

    def dfs_traverse(node, depth=0):
        if node["type"] == "mol":
            if "smiles" in node:
                try:
                    mol = Chem.MolFromSmiles(node["smiles"])
                    if mol:
                        # SMARTS for Boc group
                        boc_pattern = Chem.MolFromSmarts(
                            "[#6]-[#6](-[#6])(-[#6])-[#8]-[#6](=[#8])-[#7]"
                        )
                        if mol.HasSubstructMatch(boc_pattern):
                            boc_at_depths.add(depth)
                            print(f"Found Boc group at depth {depth}")
                except:
                    pass

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # If Boc group is present at multiple depths, it's maintained throughout synthesis
    return len(boc_at_depths) >= 2
