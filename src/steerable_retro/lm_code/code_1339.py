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
    Detects the use of halogenated aromatic compounds throughout the synthesis.
    """
    halogen_count = 0

    def dfs_traverse(node):
        nonlocal halogen_count

        if node["type"] == "mol" and "smiles" in node:
            smiles = node["smiles"]
            mol = Chem.MolFromSmiles(smiles)

            if mol:
                # Check for halogenated aromatics
                fluoro_aromatic = Chem.MolFromSmarts("[F][c]")
                chloro_aromatic = Chem.MolFromSmarts("[Cl][c]")
                bromo_aromatic = Chem.MolFromSmarts("[Br][c]")

                if mol.HasSubstructMatch(fluoro_aromatic):
                    halogen_count += len(mol.GetSubstructMatches(fluoro_aromatic))
                if mol.HasSubstructMatch(chloro_aromatic):
                    halogen_count += len(mol.GetSubstructMatches(chloro_aromatic))
                if mol.HasSubstructMatch(bromo_aromatic):
                    halogen_count += len(mol.GetSubstructMatches(bromo_aromatic))

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Found {halogen_count} halogenated aromatic positions")
    return halogen_count >= 3  # At least 3 halogen atoms on aromatic rings
