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
    Detects if the final product contains multiple nitrogen heterocycles
    (pyridine, piperidine, piperazine).
    """
    heterocycle_count = 0

    def dfs_traverse(node):
        nonlocal heterocycle_count

        if node["type"] == "mol" and not node.get(
            "children"
        ):  # Final product (leaf node)
            smiles = node["smiles"]
            mol = Chem.MolFromSmiles(smiles)

            if mol:
                # Pyridine pattern
                pyridine_pattern = Chem.MolFromSmarts("n1ccccc1")
                # Piperidine pattern
                piperidine_pattern = Chem.MolFromSmarts("N1CCCCC1")
                # Piperazine pattern
                piperazine_pattern = Chem.MolFromSmarts("N1CCNCC1")

                if mol.HasSubstructMatch(pyridine_pattern):
                    heterocycle_count += mol.GetSubstructMatches(
                        pyridine_pattern
                    ).__len__()
                if mol.HasSubstructMatch(piperidine_pattern):
                    heterocycle_count += mol.GetSubstructMatches(
                        piperidine_pattern
                    ).__len__()
                if mol.HasSubstructMatch(piperazine_pattern):
                    heterocycle_count += mol.GetSubstructMatches(
                        piperazine_pattern
                    ).__len__()

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    if heterocycle_count >= 3:
        print(f"Found {heterocycle_count} nitrogen heterocycles in final product")
        return True
    return False
