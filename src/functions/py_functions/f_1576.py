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
    Detects if the synthesis route involves an enamine intermediate.
    """
    has_enamine = False

    def dfs_traverse(node):
        nonlocal has_enamine

        if node["type"] == "mol":
            # Check for enamine pattern
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                enamine_pattern = Chem.MolFromSmarts("[#6]-[#7](-[#6])-[#6]=[#6]")
                if mol.HasSubstructMatch(enamine_pattern):
                    has_enamine = True
                    print(
                        f"Detected enamine intermediate in molecule: {node['smiles']}"
                    )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return has_enamine
