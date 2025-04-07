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
    Detects if the synthesis route involves a pyridazine core that is maintained
    throughout the synthesis.
    """
    final_product_has_pyridazine = False
    intermediate_has_pyridazine = False

    def dfs_traverse(node):
        nonlocal final_product_has_pyridazine, intermediate_has_pyridazine

        if node["type"] == "mol" and "smiles" in node:
            try:
                mol = Chem.MolFromSmiles(node["smiles"])
                pyridazine_pattern = Chem.MolFromSmarts("[n]1[n][c][c][c][c]1")

                if mol and mol.HasSubstructMatch(pyridazine_pattern):
                    # Check if this is the final product
                    if not node.get("children", []):
                        final_product_has_pyridazine = True
                        print("Final product has pyridazine core")
                    else:
                        intermediate_has_pyridazine = True
                        print("Intermediate has pyridazine core")
            except Exception as e:
                print(f"Error in SMARTS matching: {e}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return final_product_has_pyridazine and intermediate_has_pyridazine
