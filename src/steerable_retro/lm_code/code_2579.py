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
    Detects if the synthesis uses halogen-containing building blocks (F, Cl, Br, I).
    """
    found_halogen_building_blocks = False

    def dfs_traverse(node):
        nonlocal found_halogen_building_blocks

        if node["type"] == "mol" and node.get("in_stock", False):
            # This is a starting material
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                halogen_pattern = Chem.MolFromSmarts("[F,Cl,Br,I]")
                if mol.HasSubstructMatch(halogen_pattern):
                    found_halogen_building_blocks = True
                    print("Found halogen-containing building block")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return found_halogen_building_blocks
