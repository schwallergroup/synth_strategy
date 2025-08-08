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
    Detects preservation of complex aromatic systems (indole, sulfonamide) throughout synthesis.
    """
    indole_depths = []
    sulfonamide_depths = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "mol" and "smiles" in node:
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Check for indole
                indole_pattern = Chem.MolFromSmarts("[c]1[c][n][c]2[c][c][c][c][c]21")
                if mol.HasSubstructMatch(indole_pattern):
                    indole_depths.append(depth)

                # Check for sulfonamide
                sulfonamide_pattern = Chem.MolFromSmarts("[S](=[O])(=[O])[c]")
                if mol.HasSubstructMatch(sulfonamide_pattern):
                    sulfonamide_depths.append(depth)

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if both motifs are preserved across multiple depths
    indole_preserved = len(set(indole_depths)) >= 3
    sulfonamide_preserved = len(set(sulfonamide_depths)) >= 3

    if indole_preserved:
        print(f"Indole preserved across depths: {sorted(set(indole_depths))}")
    if sulfonamide_preserved:
        print(f"Sulfonamide preserved across depths: {sorted(set(sulfonamide_depths))}")

    return indole_preserved and sulfonamide_preserved
