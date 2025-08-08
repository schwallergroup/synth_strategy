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
    Detects if the synthesis involves heterocyclic compounds.
    """
    contains_heterocycles = False

    def dfs_traverse(node):
        nonlocal contains_heterocycles

        if node["type"] == "mol" and "smiles" in node:
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Check for common heterocycles
                furan_pattern = Chem.MolFromSmarts("[o]1[c][c][c][c]1")
                imidazole_pattern = Chem.MolFromSmarts("[n]1[c][n][c][c]1")
                pyridine_pattern = Chem.MolFromSmarts("[n]1[c][c][c][c][c]1")

                if (
                    mol.HasSubstructMatch(furan_pattern)
                    or mol.HasSubstructMatch(imidazole_pattern)
                    or mol.HasSubstructMatch(pyridine_pattern)
                ):
                    contains_heterocycles = True
                    print(f"Found heterocycle in molecule: {node['smiles']}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return contains_heterocycles
