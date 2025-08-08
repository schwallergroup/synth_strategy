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
    This function detects retention of aryl halide (specifically aryl iodide) throughout the synthesis.
    It checks if aryl iodide is present in all molecules along the main synthetic path.
    """
    aryl_iodide_present_in_all = True

    def dfs_traverse(node):
        nonlocal aryl_iodide_present_in_all

        if node["type"] == "mol" and "smiles" in node:
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                aryl_iodide_pattern = Chem.MolFromSmarts("[c][I]")
                if not mol.HasSubstructMatch(aryl_iodide_pattern):
                    # If this is a main product (not a reagent), check for aryl iodide
                    if not node.get("in_stock", False):
                        aryl_iodide_present_in_all = False
                        print(f"Aryl iodide not found in molecule: {node['smiles']}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return aryl_iodide_present_in_all
