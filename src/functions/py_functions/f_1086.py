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
    This function detects if the synthetic route uses an alkyne as a linker
    between aromatic rings in the final product.
    """
    has_alkyne_linker = False

    def dfs_traverse(node, depth=0):
        nonlocal has_alkyne_linker

        if node["type"] == "mol" and depth == 0:  # Final product
            try:
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol:
                    # SMARTS for aromatic-alkyne-aromatic: [c]-[C]#[C]-[c]
                    alkyne_linker_pattern = Chem.MolFromSmarts("[c]-[C]#[C]-[c]")
                    if mol.HasSubstructMatch(alkyne_linker_pattern):
                        has_alkyne_linker = True
                        print(
                            f"Detected alkyne linker in final product: {node['smiles']}"
                        )
            except:
                pass

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return has_alkyne_linker
