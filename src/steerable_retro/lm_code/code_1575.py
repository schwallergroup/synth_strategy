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
    Detects if the synthesis route involves a benzyl ether linker between aromatic rings.
    """
    has_benzyl_ether = False

    def dfs_traverse(node):
        nonlocal has_benzyl_ether

        if node["type"] == "mol":
            # Check for benzyl ether pattern
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                benzyl_ether_pattern = Chem.MolFromSmarts("c-[CH2]-[O]-c")
                if mol.HasSubstructMatch(benzyl_ether_pattern):
                    has_benzyl_ether = True
                    print(f"Detected benzyl ether linker in molecule: {node['smiles']}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return has_benzyl_ether
