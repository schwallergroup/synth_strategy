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
    This function detects if the synthesis involves building around a piperazine scaffold.
    """
    has_piperazine = False
    has_substituted_piperazine = False

    def dfs_traverse(node, depth=0):
        nonlocal has_piperazine, has_substituted_piperazine

        if node["type"] == "mol":
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Basic piperazine pattern
                piperazine_pattern = Chem.MolFromSmarts("[N]1CCN([H,C,c])CC1")

                # Disubstituted piperazine pattern
                disubst_piperazine_pattern = Chem.MolFromSmarts("[N]1CCN([C,c])CC1")

                if mol.HasSubstructMatch(piperazine_pattern):
                    has_piperazine = True

                if mol.HasSubstructMatch(disubst_piperazine_pattern):
                    has_substituted_piperazine = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Return True if we find both basic piperazine and substituted piperazine
    # indicating the synthesis involves building around this scaffold
    return has_piperazine and has_substituted_piperazine
