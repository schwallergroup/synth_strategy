#!/bin/python

"""LM-defined function for strategy description."""

import copy
import re
from collections import deque

import rdkit
import rdkit.Chem as Chem
from rdkit import Chem
from rdkit.Chem import (
    AllChem,
    Descriptors,
    Lipinski,
    rdChemReactions,
    rdFMCS,
    rdMolDescriptors,
    rdmolops,
)
from rdkit.Chem.Scaffolds import MurckoScaffold


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
