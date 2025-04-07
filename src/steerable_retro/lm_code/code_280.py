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
    This function detects if the synthetic route involves a hydrazide intermediate.
    """
    hydrazide_intermediate_found = False

    def dfs_traverse(node):
        nonlocal hydrazide_intermediate_found

        if node["type"] == "mol":
            try:
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol:
                    # Check for hydrazide group
                    hydrazide_pattern = Chem.MolFromSmarts("C(=O)NN")
                    if mol.HasSubstructMatch(hydrazide_pattern):
                        print("Hydrazide intermediate detected")
                        hydrazide_intermediate_found = True
            except:
                pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return hydrazide_intermediate_found
