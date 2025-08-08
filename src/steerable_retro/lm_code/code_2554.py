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
    Detects if the synthesis route maintains a complex polycyclic scaffold throughout.
    A complex scaffold is defined as having at least 3 rings.
    """

    def count_rings(mol):
        return mol.GetSSSR()

    scaffold_maintained = True
    min_ring_count = 3  # Threshold for "complex" polycyclic system

    def dfs_traverse(node):
        nonlocal scaffold_maintained

        if node["type"] == "mol" and node.get("smiles"):
            try:
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol:
                    ring_count = len(count_rings(mol))
                    print(f"Molecule has {ring_count} rings")
                    if ring_count < min_ring_count:
                        scaffold_maintained = False
            except:
                print("Error processing molecule SMILES")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return scaffold_maintained
