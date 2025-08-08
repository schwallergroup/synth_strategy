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
    Detects if the synthetic route uses building blocks containing heterocycles
    (specifically furan and morpholine rings)
    """
    has_furan = False
    has_morpholine = False

    def dfs_traverse(node, depth=0):
        nonlocal has_furan, has_morpholine

        if node["type"] == "mol" and "smiles" in node:
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Check for furan pattern
                furan_pattern = Chem.MolFromSmarts("o1cncc1")
                if mol.HasSubstructMatch(furan_pattern):
                    has_furan = True
                    print(f"Found furan-containing fragment at depth {depth}")

                # Check for morpholine pattern
                morpholine_pattern = Chem.MolFromSmarts("N1CCOCC1")
                if mol.HasSubstructMatch(morpholine_pattern):
                    has_morpholine = True
                    print(f"Found morpholine-containing fragment at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return has_furan and has_morpholine
