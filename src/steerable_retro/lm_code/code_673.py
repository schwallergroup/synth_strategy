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
    Detects persistence of ether linkage throughout synthesis.
    """
    # Track depths where ether linkage is present
    ether_depths = set()
    total_depths = set()

    def dfs_traverse(node, depth=0):
        nonlocal ether_depths, total_depths

        if node["type"] == "mol" and "smiles" in node:
            total_depths.add(depth)
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Check for ether linkage pattern
                pattern = Chem.MolFromSmarts("[CH2][O][CH2]")
                if mol.HasSubstructMatch(pattern):
                    ether_depths.add(depth)
                    print(f"Found ether linkage at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if ether linkage persists throughout most of the synthesis
    # (present in at least 75% of the depths)
    if total_depths:
        persistence_ratio = len(ether_depths) / len(total_depths)
        strategy_present = persistence_ratio >= 0.75
        print(f"Ether linkage persistence ratio: {persistence_ratio}")
    else:
        strategy_present = False

    return strategy_present
