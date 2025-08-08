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
    This function detects the transformation sequence: Ketone → Enaminone → Pyrazole.
    """
    # Track the presence of each intermediate in the sequence
    ketone_found = False
    enaminone_found = False
    pyrazole_found = False

    def dfs_traverse(node):
        nonlocal ketone_found, enaminone_found, pyrazole_found

        if node["type"] == "mol" and "smiles" in node:
            mol = Chem.MolFromSmiles(node["smiles"])

            if mol:
                # Check for ketone pattern
                ketone_pattern = Chem.MolFromSmarts("[C][C](=O)[c,C]")
                if mol.HasSubstructMatch(ketone_pattern):
                    ketone_found = True
                    print("Detected ketone intermediate")

                # Check for enaminone pattern
                enaminone_pattern = Chem.MolFromSmarts("CN(C)C=CC(=O)[c,C]")
                if mol.HasSubstructMatch(enaminone_pattern):
                    enaminone_found = True
                    print("Detected enaminone intermediate")

                # Check for pyrazole pattern
                pyrazole_pattern = Chem.MolFromSmarts("c1[nH]ncc1")
                if mol.HasSubstructMatch(pyrazole_pattern):
                    pyrazole_found = True
                    print("Detected pyrazole in product")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # Return True only if all three intermediates were found in the correct sequence
    return ketone_found and enaminone_found and pyrazole_found
