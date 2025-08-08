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
    Detects if the synthesis involves a tosylated pyrazole fragment.
    """
    found_tosylated_pyrazole = False

    def dfs_traverse(node):
        nonlocal found_tosylated_pyrazole

        if node["type"] == "mol" and "smiles" in node:
            try:
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol:
                    # Pattern for tosyl group
                    tosyl_pattern = Chem.MolFromSmarts(
                        "[#16](=[#8])(=[#8])-[#6]1:[#6]:[#6]:[#6](-[#6]):[#6]:[#6]:1"
                    )
                    # Pattern for pyrazole
                    pyrazole_pattern = Chem.MolFromSmarts("[#7]1:[#6]:[#7]:[#6]:[#6]:1")

                    if mol.HasSubstructMatch(tosyl_pattern) and mol.HasSubstructMatch(
                        pyrazole_pattern
                    ):
                        found_tosylated_pyrazole = True
                        print("Found tosylated pyrazole fragment")
            except:
                pass

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return found_tosylated_pyrazole
