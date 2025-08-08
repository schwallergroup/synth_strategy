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
    This function detects a synthesis strategy involving multiple nitrogen-containing
    heterocycles (pyridine, pyrimidine, piperazine, morpholine).
    """
    # Define SMARTS patterns for heterocycles
    pyridine_pattern = Chem.MolFromSmarts("c1cccnc1")
    pyrimidine_pattern = Chem.MolFromSmarts("c1cncnc1")
    piperazine_pattern = Chem.MolFromSmarts("[#7]1[#6][#6][#7][#6][#6]1")
    morpholine_pattern = Chem.MolFromSmarts("[#7]1[#6][#6][#8][#6][#6]1")

    # Track if we found the heterocycles
    found_pyridine = False
    found_pyrimidine = False
    found_piperazine = False
    found_morpholine = False

    def dfs_traverse(node):
        nonlocal found_pyridine, found_pyrimidine, found_piperazine, found_morpholine

        if node["type"] == "mol":
            try:
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol:
                    if mol.HasSubstructMatch(pyridine_pattern):
                        found_pyridine = True
                    if mol.HasSubstructMatch(pyrimidine_pattern):
                        found_pyrimidine = True
                    if mol.HasSubstructMatch(piperazine_pattern):
                        found_piperazine = True
                    if mol.HasSubstructMatch(morpholine_pattern):
                        found_morpholine = True
            except:
                pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Count how many heterocycle types we found
    heterocycle_count = sum([found_pyridine, found_pyrimidine, found_piperazine, found_morpholine])

    # Return True if at least 3 different heterocycle types were found
    if heterocycle_count >= 3:
        print(f"Found {heterocycle_count} different nitrogen heterocycle types")
        return True
    return False
