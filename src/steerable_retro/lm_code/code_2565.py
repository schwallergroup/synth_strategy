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
    This function detects if the synthesis includes installation of an ethylene linker
    connecting two heterocyclic or aromatic systems.
    """
    found_ethylene_linker = False

    def dfs_traverse(node):
        nonlocal found_ethylene_linker

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            product = rsmi.split(">")[-1]

            # Pattern for ethylene linker between aromatic/heterocyclic systems
            ethylene_linker_pattern = Chem.MolFromSmarts("c[CH2][CH2]n")

            product_mol = Chem.MolFromSmiles(product) if product else None

            if product_mol and product_mol.HasSubstructMatch(ethylene_linker_pattern):
                found_ethylene_linker = True
                print("Found ethylene linker installation")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return found_ethylene_linker
