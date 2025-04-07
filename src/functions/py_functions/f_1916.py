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
    Detects if the synthesis involves installation of a long alkyl chain (C6 or longer)
    as a linker between functional groups.
    """
    has_long_alkyl_linker = False

    def dfs_traverse(node):
        nonlocal has_long_alkyl_linker

        if node["type"] == "mol":
            smiles = node["smiles"]
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                # Look for long alkyl chains (6 or more carbons in a row)
                if mol.HasSubstructMatch(
                    Chem.MolFromSmarts("[#6]~[#6]~[#6]~[#6]~[#6]~[#6]")
                ):
                    # Check if the chain connects functional groups
                    if mol.HasSubstructMatch(
                        Chem.MolFromSmarts("[#6]~[#6]~[#6]~[#6]~[#6]~[#6][O,N,S]")
                    ) and mol.HasSubstructMatch(
                        Chem.MolFromSmarts("[O,N,S][#6]~[#6]~[#6]~[#6]~[#6]~[#6]")
                    ):
                        has_long_alkyl_linker = True
                        print("Detected long alkyl chain as linker")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return has_long_alkyl_linker
