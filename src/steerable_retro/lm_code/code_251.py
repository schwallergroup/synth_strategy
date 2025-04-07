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
    This function detects a sequential functional group interconversion strategy:
    nitro → amine → sulfonamide → N-alkylated sulfonamide
    """
    # Track if we've seen each transformation in the correct order
    seen_nitro = False
    seen_amine = False
    seen_sulfonamide = False
    seen_n_alkylated_sulfonamide = False

    # Track the depth at which each functional group appears
    nitro_depth = -1
    amine_depth = -1
    sulfonamide_depth = -1
    n_alkylated_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal seen_nitro, seen_amine, seen_sulfonamide, seen_n_alkylated_sulfonamide
        nonlocal nitro_depth, amine_depth, sulfonamide_depth, n_alkylated_depth

        if node["type"] == "mol":
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Check for nitro group
                if mol.HasSubstructMatch(Chem.MolFromSmarts("[N+](=O)[O-]")):
                    seen_nitro = True
                    nitro_depth = depth
                    print(f"Found nitro group at depth {depth}")

                # Check for primary amine
                if mol.HasSubstructMatch(Chem.MolFromSmarts("[NH2]")):
                    seen_amine = True
                    amine_depth = depth
                    print(f"Found primary amine at depth {depth}")

                # Check for sulfonamide
                if mol.HasSubstructMatch(Chem.MolFromSmarts("[NH]S(=O)(=O)[#6]")):
                    seen_sulfonamide = True
                    sulfonamide_depth = depth
                    print(f"Found sulfonamide at depth {depth}")

                # Check for N-alkylated sulfonamide
                if mol.HasSubstructMatch(
                    Chem.MolFromSmarts("[#6][N]S(=O)(=O)[#6]")
                ) and not mol.HasSubstructMatch(Chem.MolFromSmarts("[NH]S(=O)(=O)[#6]")):
                    seen_n_alkylated_sulfonamide = True
                    n_alkylated_depth = depth
                    print(f"Found N-alkylated sulfonamide at depth {depth}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if all transformations were found in the correct order
    correct_sequence = (
        seen_nitro
        and seen_amine
        and seen_sulfonamide
        and seen_n_alkylated_sulfonamide
        and nitro_depth > amine_depth > sulfonamide_depth > n_alkylated_depth
    )

    print(f"Nitro to amine to sulfonamide strategy detected: {correct_sequence}")
    return correct_sequence
