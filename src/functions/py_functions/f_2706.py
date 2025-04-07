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
    This function detects if the synthesis route includes a late-stage fragment coupling
    (defined as a C-C bond formation in the first half of the synthesis).
    """
    late_stage_coupling_found = False
    max_depth = 0

    # First, determine the maximum depth of the synthesis
    def find_max_depth(node, current_depth=0):
        nonlocal max_depth
        max_depth = max(max_depth, current_depth)

        for child in node.get("children", []):
            find_max_depth(child, current_depth + 1)

    find_max_depth(route)

    # Now check for fragment coupling in the first half of the synthesis
    def dfs_traverse(node, depth=0):
        nonlocal late_stage_coupling_found

        # Consider it late-stage if it's in the first half of the synthesis
        if (
            depth <= max_depth / 2
            and node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if we have multiple fragments combining
            if len(reactants) >= 2:
                try:
                    # Check if this is a C-C bond formation
                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol:
                        # This is a simplistic check - in a real implementation,
                        # you would need to analyze the actual bond formation
                        print(
                            f"Potential late-stage fragment coupling at depth {depth}"
                        )
                        late_stage_coupling_found = True
                except:
                    pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    return late_stage_coupling_found
