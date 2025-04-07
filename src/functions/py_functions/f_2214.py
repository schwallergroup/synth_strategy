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
    This function checks if a fluorinated aromatic ring is maintained throughout the synthesis.
    """
    all_mols_have_f_aromatic = True

    def dfs_traverse(node):
        nonlocal all_mols_have_f_aromatic

        if node["type"] == "mol" and node["smiles"]:
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol is not None:
                # Check if molecule has fluorinated aromatic ring
                has_f_aromatic = mol.HasSubstructMatch(Chem.MolFromSmarts("[c][F]"))
                if not has_f_aromatic and not node.get("in_stock", False):
                    # If it's not a starting material and doesn't have F-aromatic, mark as False
                    all_mols_have_f_aromatic = False
                    print(
                        f"Found molecule without fluorinated aromatic: {node['smiles']}"
                    )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Fluorinated aromatic maintained throughout: {all_mols_have_f_aromatic}")
    return all_mols_have_f_aromatic
