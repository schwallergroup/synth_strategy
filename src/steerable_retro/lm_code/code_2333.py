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
    This function detects if the synthesis involves heterocyclic compounds
    like benzothiazole and indole.
    """
    has_benzothiazole = False
    has_indole = False

    def dfs_traverse(node):
        nonlocal has_benzothiazole, has_indole

        if node["type"] == "mol" and "smiles" in node:
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Check for benzothiazole
                if mol.HasSubstructMatch(Chem.MolFromSmarts("c1ccc2scnc2c1")):
                    has_benzothiazole = True
                    print("Benzothiazole detected")

                # Check for indole
                if mol.HasSubstructMatch(Chem.MolFromSmarts("c1ccc2[nH]ccc2c1")):
                    has_indole = True
                    print("Indole detected")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return has_benzothiazole and has_indole
