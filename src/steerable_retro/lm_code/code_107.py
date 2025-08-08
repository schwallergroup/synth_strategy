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
    This function detects if heterocyclic scaffolds (pyridine, isoxazole) are maintained throughout synthesis.
    """
    pyridine_present = False
    isoxazole_present = False

    def dfs_traverse(node):
        nonlocal pyridine_present, isoxazole_present

        if node["type"] == "mol" and "smiles" in node:
            try:
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol:
                    # Check for pyridine
                    pyridine_pattern = Chem.MolFromSmarts("c1ccncc1")
                    # Check for isoxazole
                    isoxazole_pattern = Chem.MolFromSmarts("c1conc1")

                    if mol.HasSubstructMatch(pyridine_pattern):
                        pyridine_present = True
                    if mol.HasSubstructMatch(isoxazole_pattern):
                        isoxazole_present = True
            except:
                pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from root
    dfs_traverse(route)

    # Check if both heterocycles are present
    if pyridine_present and isoxazole_present:
        print("Detected maintenance of heterocyclic scaffolds (pyridine and isoxazole)")
        return True
    return False
