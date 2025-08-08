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
    This function detects a synthetic strategy that utilizes an aromatic amine
    as a key precursor.
    """
    aromatic_amine_found = False

    def dfs_traverse(node, depth=0):
        nonlocal aromatic_amine_found

        if node["type"] == "mol" and depth > 1:  # Look for precursors (higher depth)
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Check for aromatic amine
                aromatic_amine_pattern = Chem.MolFromSmarts("[c][NH2]")
                if mol.HasSubstructMatch(aromatic_amine_pattern):
                    aromatic_amine_found = True
                    print(f"Aromatic amine precursor found at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    if aromatic_amine_found:
        print("Aromatic amine precursor strategy detected")

    return aromatic_amine_found
