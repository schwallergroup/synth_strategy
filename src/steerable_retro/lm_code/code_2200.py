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
    This function detects if protecting groups (Boc and/or tert-butyl ester)
    persist through multiple steps of the synthesis.
    """
    # Define patterns for protecting groups
    boc_pattern = Chem.MolFromSmarts("CC(C)(C)OC(=O)N")
    tbutyl_ester_pattern = Chem.MolFromSmarts("CC(C)(C)OC(=O)C")

    # Track depths where protecting groups are present
    boc_depths = set()
    tbutyl_ester_depths = set()
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal max_depth

        max_depth = max(max_depth, depth)

        if node["type"] == "mol":
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                if mol.HasSubstructMatch(boc_pattern):
                    boc_depths.add(depth)
                if mol.HasSubstructMatch(tbutyl_ester_pattern):
                    tbutyl_ester_depths.add(depth)

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if protecting groups persist through multiple steps
    boc_persistence = len(boc_depths) >= 3  # Present in at least 3 depths
    tbutyl_persistence = len(tbutyl_ester_depths) >= 3  # Present in at least 3 depths

    print(f"Boc group present at depths: {sorted(boc_depths)}")
    print(f"tert-Butyl ester present at depths: {sorted(tbutyl_ester_depths)}")
    print(f"Boc persistence: {boc_persistence}")
    print(f"tert-Butyl ester persistence: {tbutyl_persistence}")

    return boc_persistence or tbutyl_persistence
