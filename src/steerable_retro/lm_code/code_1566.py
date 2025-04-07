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
    This function detects if the synthesis involves construction of a polycyclic scaffold
    (molecule with multiple rings in the final product).
    """
    polycyclic_found = False

    def dfs_traverse(node, depth=0):
        nonlocal polycyclic_found

        if node["type"] == "mol" and depth == 0:  # Final product
            if "smiles" in node:
                smiles = node["smiles"]
                mol = Chem.MolFromSmiles(smiles)
                if mol is not None:
                    ring_count = Chem.rdMolDescriptors.CalcNumRings(mol)
                    if ring_count >= 3:
                        print(f"Polycyclic scaffold detected with {ring_count} rings")
                        polycyclic_found = True

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    return polycyclic_found
