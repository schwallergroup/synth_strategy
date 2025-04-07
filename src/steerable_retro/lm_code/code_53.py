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
    This function detects a strategy involving both indole and piperazine scaffolds
    in the final product.
    """
    has_indole = False
    has_piperazine = False

    def dfs_traverse(node):
        nonlocal has_indole, has_piperazine

        if node["type"] == "mol" and node.get("smiles"):
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Check for indole scaffold
                indole_pattern = Chem.MolFromSmarts("c1ccc2[nH]ccc2c1")
                if mol.HasSubstructMatch(indole_pattern):
                    has_indole = True
                    print("Indole scaffold detected")

                # Check for piperazine scaffold
                piperazine_pattern = Chem.MolFromSmarts("N1CCN([#6])CC1")
                if mol.HasSubstructMatch(piperazine_pattern):
                    has_piperazine = True
                    print("Piperazine scaffold detected")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if both scaffolds are present
    return has_indole and has_piperazine
