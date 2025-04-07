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
    This function detects a strategy using ester intermediates that are later
    converted to other functional groups.
    """
    ester_depths = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "mol" and "smiles" in node:
            try:
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol:
                    # Check for ester group
                    ester_pattern = Chem.MolFromSmarts("[#6]-[#8]-[#6](=[#8])-[#6]")
                    if mol.HasSubstructMatch(ester_pattern):
                        ester_depths.append(depth)
                        print(f"Ester group found at depth {depth}")
            except:
                print(f"Error processing molecule SMILES at depth {depth}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if esters are present in intermediate steps but not in final product
    if ester_depths and 0 not in ester_depths:
        print("Ester intermediate strategy detected")
        return True
    return False
