#!/bin/python

"""LM-defined function for strategy description."""

from rdkit.Chem import AllChem, rdFMCS
import copy
from collections import deque
import rdkit.Chem as Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdChemReactions
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
import rdkit.Chem.rdFMCS
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem import AllChem, Descriptors, Lipinski
from rdkit.Chem import rdmolops
import re
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import AllChem, Descriptors
import traceback
import rdkit
from collections import Counter


def main(route):
    """
    This function detects a synthetic strategy where the number of rings progressively
    increases throughout the synthesis.
    """
    ring_counts_by_depth = {}

    def dfs_traverse(node, depth=0):
        if node["type"] == "mol" and "smiles" in node:
            try:
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol:
                    ring_info = mol.GetRingInfo()
                    ring_count = ring_info.NumRings()

                    if depth not in ring_counts_by_depth:
                        ring_counts_by_depth[depth] = ring_count
                    else:
                        ring_counts_by_depth[depth] = max(ring_counts_by_depth[depth], ring_count)
            except:
                print(f"Error processing molecule SMILES: {node['smiles']}")

        # Process children with increased depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if ring count increases as depth increases
    depths = sorted(ring_counts_by_depth.keys())
    if len(depths) < 2:
        return False

    is_increasing = True
    for i in range(1, len(depths)):
        if ring_counts_by_depth[depths[i - 1]] <= ring_counts_by_depth[depths[i]]:
            is_increasing = False
            break

    if is_increasing:
        print(f"Detected progressive ring count increase: {ring_counts_by_depth}")

    return is_increasing
