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
    This function detects if the synthetic route involves heterocyclic structures
    like benzimidazole and pyridine.
    """
    heterocycles_found = False

    def dfs_traverse(node):
        nonlocal heterocycles_found

        if node["type"] == "mol" and "smiles" in node:
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Check for benzimidazole
                benzimidazole_pattern = Chem.MolFromSmarts("c1nc2ccccc2[nH]1")
                # Check for pyridine
                pyridine_pattern = Chem.MolFromSmarts("c1ccncc1")

                if mol.HasSubstructMatch(benzimidazole_pattern) or mol.HasSubstructMatch(
                    pyridine_pattern
                ):
                    print("Heterocyclic structure detected")
                    heterocycles_found = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return heterocycles_found
