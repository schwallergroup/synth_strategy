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
    This function detects a strategy where an α,β-unsaturated carbonyl system
    is preserved throughout the synthesis.
    """
    # SMARTS pattern for α,β-unsaturated carbonyl
    unsaturated_carbonyl_pattern = Chem.MolFromSmarts("[C]=[C]C(=O)")

    # Track if all non-starting materials have the unsaturated carbonyl
    all_intermediates_have_pattern = True
    found_any_intermediate = False

    def dfs_traverse(node):
        nonlocal all_intermediates_have_pattern, found_any_intermediate

        if node["type"] == "mol" and not node.get("in_stock", False):
            found_any_intermediate = True
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                if not mol.HasSubstructMatch(unsaturated_carbonyl_pattern):
                    all_intermediates_have_pattern = False
                    print(f"Found intermediate without α,β-unsaturated carbonyl: {node['smiles']}")
                else:
                    print(f"Found intermediate with α,β-unsaturated carbonyl: {node['smiles']}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    strategy_present = all_intermediates_have_pattern and found_any_intermediate
    print(f"Preserved α,β-unsaturated carbonyl strategy detected: {strategy_present}")
    return strategy_present
