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
    Detects if a nitro-aromatic group is maintained throughout the synthesis.
    """
    nitro_aromatic_depths = set()

    def dfs_traverse(node, depth=0):
        if node["type"] == "mol" and "smiles" in node:
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Check for nitro-aromatic group
                nitro_aromatic_pattern = Chem.MolFromSmarts("c[N+](=[O])[O-]")
                if mol.HasSubstructMatch(nitro_aromatic_pattern):
                    nitro_aromatic_depths.add(depth)

        # Continue traversing the tree
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if nitro-aromatic group is present at multiple depths
    strategy_present = len(nitro_aromatic_depths) >= 3  # Present at multiple steps
    print(f"Nitro-aromatic maintained strategy detected: {strategy_present}")
    print(f"Found at depths: {nitro_aromatic_depths}")
    return strategy_present
