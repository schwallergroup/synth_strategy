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
    Detects synthesis routes that utilize a dihalogenated aryl fragment (containing both Cl and Br).
    """
    dihalogenated_fragment_used = False

    def dfs_traverse(node):
        nonlocal dihalogenated_fragment_used

        if node["type"] == "mol":
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Check for both Cl and Br on aromatic rings
                cl_pattern = Chem.MolFromSmarts("[c][Cl]")
                br_pattern = Chem.MolFromSmarts("[c][Br]")

                if mol.HasSubstructMatch(cl_pattern) and mol.HasSubstructMatch(br_pattern):
                    dihalogenated_fragment_used = True
                    print("Detected dihalogenated aryl fragment (Cl and Br)")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return dihalogenated_fragment_used
