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
    This function detects if the synthesis maintains a trifluoromethyl group
    throughout the synthesis.
    """
    steps_with_cf3 = 0
    total_steps = 0

    def dfs_traverse(node, depth=0):
        nonlocal steps_with_cf3, total_steps

        if node["type"] == "reaction":
            total_steps += 1
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                product = rsmi.split(">")[-1]

                # Check for trifluoromethyl pattern
                cf3_pattern = Chem.MolFromSmarts("[c]-[C](-[F])(-[F])-[F]")

                product_mol = Chem.MolFromSmiles(product) if product else None

                if product_mol and product_mol.HasSubstructMatch(cf3_pattern):
                    print("Found trifluoromethyl group at depth", depth)
                    steps_with_cf3 += 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)

    # Return True if CF3 is present in all steps
    return steps_with_cf3 > 0 and steps_with_cf3 == total_steps
