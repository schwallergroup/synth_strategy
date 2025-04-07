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


def main(route, threshold=3):
    """
    Detects if the final product in the synthesis route has multiple aromatic rings
    """
    aromatic_ring_count = 0

    def dfs_traverse(node, depth=0):
        nonlocal aromatic_ring_count

        if node["type"] == "mol" and depth == 0:  # Final product
            print(f"Analyzing final product: {node['smiles']}")
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Count aromatic rings
                aromatic_rings = 0
                ring_info = mol.GetRingInfo()
                for ring in ring_info.AtomRings():
                    # Check if all atoms in the ring are aromatic
                    if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
                        aromatic_rings += 1
                        print(f"Found aromatic ring: {[idx for idx in ring]}")

                aromatic_ring_count = aromatic_rings
                print(f"Found {aromatic_ring_count} aromatic rings in final product")
            else:
                print(f"Failed to parse molecule SMILES: {node['smiles']}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return aromatic_ring_count >= threshold
