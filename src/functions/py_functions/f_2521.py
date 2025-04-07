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
    This function detects if the fluorobenzene aromatic scaffold is preserved throughout the synthesis.
    """
    all_intermediates_have_scaffold = True

    def dfs_traverse(node):
        nonlocal all_intermediates_have_scaffold

        if node["type"] == "mol" and not node.get("in_stock", False):
            if node.get("smiles"):
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol:
                    # Check for fluorobenzene scaffold
                    fluoro_benzene = Chem.MolFromSmarts("c1c(F)cccc1")
                    if not mol.HasSubstructMatch(fluoro_benzene):
                        all_intermediates_have_scaffold = False
                        print(
                            "Molecule without fluorobenzene scaffold:", node["smiles"]
                        )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    print(f"Aromatic scaffold preservation detected: {all_intermediates_have_scaffold}")
    return all_intermediates_have_scaffold
