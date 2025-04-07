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
    This function detects if a synthetic route maintains an α,β-unsaturated ester motif
    throughout multiple steps of the synthesis.
    """
    steps_with_motif = 0
    total_steps = 0

    def dfs_traverse(node):
        nonlocal steps_with_motif, total_steps

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            total_steps += 1
            rsmi = node["metadata"]["rsmi"]
            product = rsmi.split(">")[-1]

            # Check for α,β-unsaturated ester in product
            unsaturated_ester_pattern = Chem.MolFromSmarts("[C]=[C][C](=[O])[O][C]")

            if product.strip():
                try:
                    mol = Chem.MolFromSmiles(product)
                    if mol and mol.HasSubstructMatch(unsaturated_ester_pattern):
                        steps_with_motif += 1
                        print(f"Found α,β-unsaturated ester in step: {rsmi}")
                except:
                    pass

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if the motif is present in most steps (>= 75%)
    result = total_steps > 0 and steps_with_motif / total_steps >= 0.75
    print(
        f"α,β-unsaturated ester strategy detected: {result} (present in {steps_with_motif}/{total_steps} steps)"
    )
    return result
