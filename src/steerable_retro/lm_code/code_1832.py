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
    This function detects a strategy involving a long aliphatic chain (C12+)
    with carboxylic acid functionalities.
    """
    has_long_chain = False
    has_carboxylic_acid = False

    def dfs_traverse(node, depth=0):
        nonlocal has_long_chain, has_carboxylic_acid

        if node["type"] == "mol":
            mol = Chem.MolFromSmiles(node["smiles"])
            if not mol:
                return

            # Check for carboxylic acid
            if mol.HasSubstructMatch(Chem.MolFromSmarts("C(=O)[OH]")):
                has_carboxylic_acid = True
                print(f"Found carboxylic acid at depth {depth}")

            # Check for long aliphatic chain (C12+)
            # This is a simplified approach - in practice you might need a more sophisticated algorithm
            aliphatic_chain_pattern = Chem.MolFromSmarts("CCCCCCCCCCCCC")  # C13 chain
            if mol.HasSubstructMatch(aliphatic_chain_pattern):
                has_long_chain = True
                print(f"Found long aliphatic chain at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    strategy_present = has_long_chain and has_carboxylic_acid
    print(f"Long chain carboxylic acid strategy detected: {strategy_present}")
    return strategy_present
