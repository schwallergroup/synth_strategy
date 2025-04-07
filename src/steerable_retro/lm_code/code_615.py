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
    This function detects a strategy involving incorporation of a trifluoromethoxy group.
    """
    found_trifluoromethoxy = False

    def dfs_traverse(node, depth=0):
        nonlocal found_trifluoromethoxy

        if node["type"] == "mol":
            # Check if molecule contains trifluoromethoxy group
            smiles = node["smiles"]
            mol = Chem.MolFromSmiles(smiles)

            if mol:
                # SMARTS pattern for trifluoromethoxy group
                trifluoromethoxy_pattern = Chem.MolFromSmarts("[O][C]([F])([F])[F]")

                if mol.HasSubstructMatch(trifluoromethoxy_pattern):
                    found_trifluoromethoxy = True
                    print(f"Found trifluoromethoxy group at depth {depth}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return found_trifluoromethoxy
