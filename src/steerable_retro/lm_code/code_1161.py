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
    This function detects if a benzimidazole ring is formed in a late stage of the synthesis.
    Late stage is defined as the first or second reaction from the target molecule.
    """
    # Define the benzimidazole pattern
    benzimidazole_pattern = Chem.MolFromSmarts("c1nc2ccccc2n1")

    # Track at which depth the benzimidazole is formed
    benzimidazole_formation_depth = None
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal benzimidazole_formation_depth, max_depth

        max_depth = max(max_depth, depth)

        if node["type"] == "mol":
            # Check if this molecule contains a benzimidazole
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol and mol.HasSubstructMatch(benzimidazole_pattern):
                # If we haven't seen a benzimidazole formation yet, or this is at a lower depth
                if benzimidazole_formation_depth is None or depth < benzimidazole_formation_depth:
                    benzimidazole_formation_depth = depth

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if benzimidazole is formed in a late stage (first or second reaction)
    is_late_stage = benzimidazole_formation_depth is not None and benzimidazole_formation_depth <= 1

    print(f"Benzimidazole formation detected at depth: {benzimidazole_formation_depth}")
    print(f"Maximum depth of synthesis: {max_depth}")
    print(f"Is late-stage benzimidazole formation: {is_late_stage}")

    return is_late_stage
