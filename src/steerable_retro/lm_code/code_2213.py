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


def main(route, min_rings=3):
    """
    Detects if the final product contains at least the specified number of aromatic rings.
    """
    ring_count = 0

    def count_aromatic_rings(mol):
        if not mol:
            return 0

        # Find all aromatic rings
        ring_info = mol.GetRingInfo().AtomRings()
        aromatic_ring_count = 0

        for ring in ring_info:
            is_aromatic = True
            for atom_idx in ring:
                atom = mol.GetAtomWithIdx(atom_idx)
                if not atom.GetIsAromatic():
                    is_aromatic = False
                    break
            if is_aromatic:
                aromatic_ring_count += 1
                print(f"Found aromatic ring: {[i for i in ring]}")

        return aromatic_ring_count

    def dfs_traverse(node, depth=0):
        nonlocal ring_count

        if depth == 0 and node["type"] == "mol" and "smiles" in node:
            print(f"Analyzing final product with SMILES: {node['smiles']}")
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                ring_count = count_aromatic_rings(mol)
                print(f"Found {ring_count} aromatic rings in final product")
            else:
                print(f"Error: Could not parse SMILES string: {node['smiles']}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    result = ring_count >= min_rings
    print(f"Result: {result} (Found {ring_count} rings, minimum required: {min_rings})")
    return result
