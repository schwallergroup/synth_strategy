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
    Detects if the synthesis involves a heterocycle (like thiophene) that is maintained
    throughout the synthesis.
    """
    has_heterocycle = False

    def dfs_traverse(node):
        nonlocal has_heterocycle

        if node["type"] == "mol" and "smiles" in node:
            smiles = node["smiles"]
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                # Check for thiophene pattern
                thiophene_pattern = Chem.MolFromSmarts("c1cscc1")
                if mol.HasSubstructMatch(thiophene_pattern):
                    has_heterocycle = True
                    print("Found thiophene heterocycle")

                # Check for other common heterocycles if needed
                # furan, pyrrole, etc.

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(f"Heterocycle-containing synthesis detected: {has_heterocycle}")
    return has_heterocycle
