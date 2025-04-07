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
    Detects if an iodine atom is preserved throughout the entire synthesis.
    """
    all_mols_have_iodine = True
    mol_count = 0

    def dfs_traverse(node):
        nonlocal all_mols_have_iodine, mol_count

        if node["type"] == "mol" and not node.get("in_stock", False):
            if node["smiles"]:
                mol_count += 1
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol:
                    iodine_pattern = Chem.MolFromSmarts("[#53]")
                    if not mol.HasSubstructMatch(iodine_pattern):
                        all_mols_have_iodine = False

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    result = all_mols_have_iodine and mol_count > 0
    print(f"Iodine preservation throughout synthesis: {result}")
    return result
