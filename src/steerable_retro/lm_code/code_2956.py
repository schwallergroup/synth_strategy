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
    This function detects if the synthetic route uses cyclic amine building blocks like pyrrolidine and cyclobutylamine.
    """
    pyrrolidine_found = False
    cyclobutylamine_found = False

    def dfs_traverse(node):
        nonlocal pyrrolidine_found, cyclobutylamine_found

        if node["type"] == "mol" and node.get("in_stock", False):
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Pattern for pyrrolidine
                pyrrolidine_pattern = Chem.MolFromSmarts("[N]1[C][C][C][C]1")
                # Pattern for cyclobutylamine
                cyclobutylamine_pattern = Chem.MolFromSmarts("[N][C]1[C][C][C]1")

                if mol.HasSubstructMatch(pyrrolidine_pattern):
                    pyrrolidine_found = True
                    print(f"Found pyrrolidine building block: {node['smiles']}")

                if mol.HasSubstructMatch(cyclobutylamine_pattern):
                    cyclobutylamine_found = True
                    print(f"Found cyclobutylamine building block: {node['smiles']}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    result = pyrrolidine_found or cyclobutylamine_found
    print(f"Cyclic amine building blocks detected: {result}")
    return result
