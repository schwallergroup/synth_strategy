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
    This function detects a synthetic strategy involving a thioether-containing linker.
    """
    found_thioether_linker = False

    def dfs_traverse(node):
        nonlocal found_thioether_linker

        if node["type"] == "mol" and "smiles" in node:
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol is not None:
                # Look for thioether pattern with carbon chains on both sides
                thioether_pattern = Chem.MolFromSmarts("[#6][S][#6][#6]")
                if len(mol.GetSubstructMatches(thioether_pattern)) > 0:
                    found_thioether_linker = True
                    print(f"Detected thioether-containing linker in molecule: {node['smiles']}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Thioether-containing linker detected: {found_thioether_linker}")
    return found_thioether_linker
