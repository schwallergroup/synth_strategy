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
    Detects the use of TMS (trimethylsilyl) protection in the synthesis.
    """
    has_tms_protection = False

    def dfs_traverse(node):
        nonlocal has_tms_protection

        if node["type"] == "mol":
            if "smiles" in node:
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol:
                    # TMS pattern
                    tms_patt = Chem.MolFromSmarts("[#6][#14]([#6])([#6])[#6]")
                    if mol.HasSubstructMatch(tms_patt):
                        print("Found TMS protected molecule")
                        has_tms_protection = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return has_tms_protection
