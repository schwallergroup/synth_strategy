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
    This function detects if a pyrimidine-triazole biaryl system is constructed
    during the synthesis.
    """
    biaryl_system_present = False

    def dfs_traverse(node):
        nonlocal biaryl_system_present

        if node["type"] == "mol":
            if "smiles" in node:
                smiles = node["smiles"]

                # Check for pyrimidine-triazole biaryl system
                biaryl_pattern = Chem.MolFromSmarts("[n]1[n][n][c](-[c]2[n][c][c][n][c]2)[n]1")
                try:
                    mol = Chem.MolFromSmiles(smiles)
                    if mol and mol.HasSubstructMatch(biaryl_pattern):
                        biaryl_system_present = True
                        print(f"Found pyrimidine-triazole biaryl system: {smiles}")
                except:
                    pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    return biaryl_system_present
