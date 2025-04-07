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
    Detects a synthesis that combines dimethoxyphenyl and fluorophenyl moieties.
    """
    # Track if both moieties are present in the final product
    dimethoxyphenyl_present = False
    fluorophenyl_present = False

    def dfs_traverse(node, depth=0):
        nonlocal dimethoxyphenyl_present, fluorophenyl_present

        if node["type"] == "mol" and depth == 0:  # Final product
            if node.get("smiles"):
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol:
                    if mol.HasSubstructMatch(
                        Chem.MolFromSmarts("[c]1([O][C])[c]([O][C])[c][c][c][c]1")
                    ):
                        dimethoxyphenyl_present = True
                        print("Found dimethoxyphenyl moiety in final product")

                    if mol.HasSubstructMatch(
                        Chem.MolFromSmarts("[c]1[c][c][c]([F])[c][c]1")
                    ):
                        fluorophenyl_present = True
                        print("Found fluorophenyl moiety in final product")

        # Traverse children with incremented depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return dimethoxyphenyl_present and fluorophenyl_present
