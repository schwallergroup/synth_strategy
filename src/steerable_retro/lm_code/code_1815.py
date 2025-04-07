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
    Detects a synthesis involving halogenated aromatic compounds
    (specifically with bromine and fluorine substituents).
    """
    has_bromo_aromatic = False
    has_fluoro_aromatic = False

    def dfs_traverse(node):
        nonlocal has_bromo_aromatic, has_fluoro_aromatic

        if node["type"] == "mol":
            # Check for halogenated aromatics
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                bromo_aromatic_pattern = Chem.MolFromSmarts("c[Br]")
                fluoro_aromatic_pattern = Chem.MolFromSmarts("c[F]")

                if mol.HasSubstructMatch(bromo_aromatic_pattern):
                    has_bromo_aromatic = True
                    print("Found brominated aromatic compound")

                if mol.HasSubstructMatch(fluoro_aromatic_pattern):
                    has_fluoro_aromatic = True
                    print("Found fluorinated aromatic compound")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return has_bromo_aromatic and has_fluoro_aromatic
