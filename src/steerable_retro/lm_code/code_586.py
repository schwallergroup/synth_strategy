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
    Detects a strategy for synthesizing a conjugated dye system with a PEG spacer.
    """
    # Track key structural elements
    has_styryl_bridge = False
    has_peg_linker = False
    has_dimethylamino = False

    def dfs_traverse(node, depth=0):
        nonlocal has_styryl_bridge, has_peg_linker, has_dimethylamino

        if node["type"] == "mol" and depth == 0:  # Final product
            mol = Chem.MolFromSmiles(node["smiles"])
            if not mol:
                return

            # Check for styryl bridge (conjugated C=C between aromatics)
            styryl_pattern = Chem.MolFromSmarts("c-C=C-c")
            if mol.HasSubstructMatch(styryl_pattern):
                has_styryl_bridge = True
                print("Found styryl bridge in final product")

            # Check for PEG linker
            peg_pattern = Chem.MolFromSmarts("O-C-C-O-C-C-O")
            if mol.HasSubstructMatch(peg_pattern):
                has_peg_linker = True
                print("Found PEG linker in final product")

            # Check for dimethylamino group
            dimethylamino_pattern = Chem.MolFromSmarts("c-N(-C)(-C)")
            if mol.HasSubstructMatch(dimethylamino_pattern):
                has_dimethylamino = True
                print("Found dimethylamino group in final product")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Return True if we found all key elements of the conjugated dye
    return has_styryl_bridge and has_peg_linker and has_dimethylamino
