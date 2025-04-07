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
    This function detects synthesis involving furan-containing heterocyclic systems.
    """
    furan_system_detected = False

    def dfs_traverse(node):
        nonlocal furan_system_detected

        if node["type"] == "mol":
            smiles = node["smiles"]
            mol = Chem.MolFromSmiles(smiles)

            if mol:
                # Pattern for furan-containing system with amide: c1nc2c(cc1)OCC(=O)N2
                furan_pattern = Chem.MolFromSmarts("c1nc2c(cc1)OCC(=O)N2")

                if mol.HasSubstructMatch(furan_pattern):
                    print(f"Furan-containing system detected: {smiles}")
                    furan_system_detected = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return furan_system_detected
