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
    Detects a synthetic strategy where a nitrile group is present
    throughout the synthesis and retained in the final product.
    """
    nitrile_in_final_product = False
    nitrile_in_intermediates = False

    def dfs_traverse(node, depth=0):
        nonlocal nitrile_in_final_product, nitrile_in_intermediates

        if node["type"] == "mol" and "smiles" in node:
            smiles = node["smiles"]
            mol = Chem.MolFromSmiles(smiles)

            if mol and mol.HasSubstructMatch(Chem.MolFromSmarts("C#N")):
                if depth == 0:  # Final product
                    nitrile_in_final_product = True
                    print("Found nitrile in final product")
                else:  # Intermediate
                    nitrile_in_intermediates = True
                    print("Found nitrile in intermediate")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if nitrile is retained throughout synthesis
    return nitrile_in_final_product and nitrile_in_intermediates
