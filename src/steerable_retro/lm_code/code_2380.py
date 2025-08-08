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
    Detects synthesis routes that involve multiple heterocyclic systems,
    specifically thiophene, pyrimidine, morpholine, and piperazine.
    """
    # Count of different heterocycles found
    heterocycle_count = 0
    found_thiophene = False
    found_pyrimidine = False
    found_morpholine = False
    found_piperazine = False

    def dfs_traverse(node):
        nonlocal heterocycle_count, found_thiophene, found_pyrimidine, found_morpholine, found_piperazine

        if node["type"] == "mol" and "smiles" in node:
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Check for each heterocycle
                thiophene = Chem.MolFromSmarts("[#6]1[#6][#6][#16][#6]1")
                pyrimidine = Chem.MolFromSmarts("[#6]1[#7][#6][#7][#6][#6]1")
                morpholine = Chem.MolFromSmarts("[#6]1[#6][#8][#6][#6][#7]1")
                piperazine = Chem.MolFromSmarts("[#6]1[#6][#7][#6][#6][#7]1")

                if not found_thiophene and mol.HasSubstructMatch(thiophene):
                    found_thiophene = True
                    heterocycle_count += 1
                    print("Found thiophene heterocycle")

                if not found_pyrimidine and mol.HasSubstructMatch(pyrimidine):
                    found_pyrimidine = True
                    heterocycle_count += 1
                    print("Found pyrimidine heterocycle")

                if not found_morpholine and mol.HasSubstructMatch(morpholine):
                    found_morpholine = True
                    heterocycle_count += 1
                    print("Found morpholine heterocycle")

                if not found_piperazine and mol.HasSubstructMatch(piperazine):
                    found_piperazine = True
                    heterocycle_count += 1
                    print("Found piperazine heterocycle")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if at least 3 different heterocycles are found
    return heterocycle_count >= 3
