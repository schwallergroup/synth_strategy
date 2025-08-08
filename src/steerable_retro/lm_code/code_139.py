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
    This function detects a protection-deprotection strategy involving
    a benzylic amine protected as a sulfonamide.
    """
    # Track if we found both amine and sulfonamide
    found_amine = False
    found_sulfonamide = False

    # SMARTS patterns
    amine_pattern = Chem.MolFromSmarts("[c][CH2][NH2]")
    sulfonamide_pattern = Chem.MolFromSmarts("[c][CH2][NH][S](=O)(=O)[c]")

    def dfs_traverse(node):
        nonlocal found_amine, found_sulfonamide

        if node["type"] == "mol":
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                if mol.HasSubstructMatch(amine_pattern):
                    found_amine = True
                if mol.HasSubstructMatch(sulfonamide_pattern):
                    found_sulfonamide = True

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Found amine: {found_amine}")
    print(f"Found sulfonamide: {found_sulfonamide}")

    # Return True if both amine and sulfonamide are found
    return found_amine and found_sulfonamide
