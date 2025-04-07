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
    Detects if the synthesis uses halogen-containing building blocks
    (both trifluoromethyl and bromo-aromatic)
    """
    # Initialize flags
    has_trifluoromethyl = False
    has_bromo_aromatic = False

    def dfs_traverse(node):
        nonlocal has_trifluoromethyl, has_bromo_aromatic

        if node["type"] == "mol" and node.get("in_stock", False):
            # This is a starting material
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                if mol.HasSubstructMatch(Chem.MolFromSmarts("[C]([F])([F])[F]")):
                    print("Found trifluoromethyl group in starting material")
                    has_trifluoromethyl = True

                if mol.HasSubstructMatch(Chem.MolFromSmarts("[c][Br]")):
                    print("Found bromo-aromatic group in starting material")
                    has_bromo_aromatic = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if both halogen-containing groups are found in starting materials
    return has_trifluoromethyl and has_bromo_aromatic
