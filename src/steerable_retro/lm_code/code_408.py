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
    Detects if the synthetic route involves multiple heterocyclic systems
    such as pyrazole, isoxazole, and thiophene.
    """
    heterocycle_count = 0
    heterocycle_types = set()

    def dfs_traverse(node, depth=0):
        nonlocal heterocycle_count, heterocycle_types

        if node["type"] == "mol" and "smiles" in node:
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Check for different heterocycles
                pyrazole_pattern = Chem.MolFromSmarts("[n]1[n]cc[c]1")
                isoxazole_pattern = Chem.MolFromSmarts("[n]1oc[c]c1")
                thiophene_pattern = Chem.MolFromSmarts("[s]1ccc[c]1")
                morpholine_pattern = Chem.MolFromSmarts("[N]1CCOC[C]1")

                if mol.HasSubstructMatch(pyrazole_pattern):
                    heterocycle_types.add("pyrazole")
                if mol.HasSubstructMatch(isoxazole_pattern):
                    heterocycle_types.add("isoxazole")
                if mol.HasSubstructMatch(thiophene_pattern):
                    heterocycle_types.add("thiophene")
                if mol.HasSubstructMatch(morpholine_pattern):
                    heterocycle_types.add("morpholine")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    heterocycle_count = len(heterocycle_types)

    print(f"Found {heterocycle_count} different heterocycle types: {', '.join(heterocycle_types)}")
    return heterocycle_count >= 3  # Route has at least 3 different heterocycle types
