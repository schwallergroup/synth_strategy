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
    Detects a strategy involving sulfonyl chloride activation through
    perfluorophenyl sulfonate to form a sulfonamide.
    """
    # Track if we find each step in the sequence
    found_sulfonyl_chloride = False
    found_perfluorophenyl_sulfonate = False
    found_sulfonamide = False

    # SMARTS patterns
    sulfonyl_chloride_pattern = Chem.MolFromSmarts("[#16](=[#8])(=[#8])[Cl]")
    perfluorophenyl_sulfonate_pattern = Chem.MolFromSmarts(
        "[#16](=[#8])(=[#8])[#8][c]1[c]([F])[c]([F])[c]([F])[c]([F])[c]1[F]"
    )
    sulfonamide_pattern = Chem.MolFromSmarts("[#16](=[#8])(=[#8])[#7]")

    def dfs_traverse(node):
        nonlocal found_sulfonyl_chloride, found_perfluorophenyl_sulfonate, found_sulfonamide

        if node["type"] == "mol":
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Check for each functional group
                if mol.HasSubstructMatch(sulfonyl_chloride_pattern):
                    found_sulfonyl_chloride = True
                if mol.HasSubstructMatch(perfluorophenyl_sulfonate_pattern):
                    found_perfluorophenyl_sulfonate = True
                if mol.HasSubstructMatch(sulfonamide_pattern):
                    found_sulfonamide = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Check if we found the complete sequence
    result = (
        found_sulfonyl_chloride
        and found_perfluorophenyl_sulfonate
        and found_sulfonamide
    )

    print(f"Sulfonyl activation strategy detected: {result}")
    print(f"  - Found sulfonyl chloride: {found_sulfonyl_chloride}")
    print(f"  - Found perfluorophenyl sulfonate: {found_perfluorophenyl_sulfonate}")
    print(f"  - Found sulfonamide: {found_sulfonamide}")

    return result
