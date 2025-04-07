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
    Detects if the synthetic route incorporates a trifluoromethyl-containing fragment.
    """
    trifluoromethyl_pattern = Chem.MolFromSmarts("[#6]-[C]([F])([F])[F]")
    trifluoromethyl_depths = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "mol" and node.get("smiles"):
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol and mol.HasSubstructMatch(trifluoromethyl_pattern):
                trifluoromethyl_depths.append(depth)
                print(f"Found trifluoromethyl group at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if trifluoromethyl group is introduced in a specific reaction
    if len(trifluoromethyl_depths) > 0:
        # Find the earliest appearance (highest depth)
        earliest_appearance = max(trifluoromethyl_depths)
        print(f"Trifluoromethyl group first appears at depth {earliest_appearance}")

        # If it appears after depth 0, it was introduced during synthesis
        if earliest_appearance > 0:
            print("Trifluoromethyl fragment was incorporated during synthesis")
            return True

    print("No trifluoromethyl fragment incorporation detected")
    return False
