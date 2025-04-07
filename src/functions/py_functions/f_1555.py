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
    Detects a synthesis strategy where ketone and hydroxyl functional groups
    are preserved throughout the synthesis.
    """
    # Track functional groups in each molecule
    all_mols = []

    # Define SMARTS patterns
    ketone_pattern = Chem.MolFromSmarts("[#6][C](=[O])[#6]")
    hydroxyl_pattern = Chem.MolFromSmarts("[#6][OH]")

    def dfs_traverse(node, depth=0):
        if node["type"] == "mol" and node.get("smiles"):
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                all_mols.append((mol, depth))

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Sort molecules by depth (from final product to starting materials)
    all_mols.sort(key=lambda x: x[1])

    # Check if ketones and hydroxyls are preserved
    has_ketone = [bool(mol.HasSubstructMatch(ketone_pattern)) for mol, _ in all_mols]
    has_hydroxyl = [
        bool(mol.HasSubstructMatch(hydroxyl_pattern)) for mol, _ in all_mols
    ]

    # Strategy is present if both functional groups appear in the final product
    # and are preserved in at least one molecule throughout the synthesis
    ketone_preserved = has_ketone[0] and any(has_ketone[1:])
    hydroxyl_preserved = has_hydroxyl[0] and any(has_hydroxyl[1:])
    strategy_present = ketone_preserved and hydroxyl_preserved

    print(f"Preserved functional groups strategy detection:")
    print(f"  Ketone present in molecules: {has_ketone}")
    print(f"  Hydroxyl present in molecules: {has_hydroxyl}")
    print(f"  Ketone preserved: {ketone_preserved}")
    print(f"  Hydroxyl preserved: {hydroxyl_preserved}")
    print(f"  Strategy present: {strategy_present}")

    return strategy_present
