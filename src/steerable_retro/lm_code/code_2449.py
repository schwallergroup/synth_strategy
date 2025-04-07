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
    This function detects if the synthesis involves a difluorobenzyl group.
    """
    has_difluorobenzyl = False

    def dfs_traverse(node):
        nonlocal has_difluorobenzyl

        if node["type"] == "mol" and "smiles" in node:
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Define difluorobenzyl pattern
                difluorobenzyl_pattern = Chem.MolFromSmarts("[c]1[c]([F])[c]([F])[c][c][c]1[C]")
                if mol.HasSubstructMatch(difluorobenzyl_pattern):
                    has_difluorobenzyl = True
                    print(f"Found difluorobenzyl group in molecule: {node['smiles']}")

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(
        f"Synthesis {'contains' if has_difluorobenzyl else 'does not contain'} difluorobenzyl group"
    )
    return has_difluorobenzyl
