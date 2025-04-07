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
    This function detects a synthesis route that maintains a trifluoromethyl group
    throughout the synthesis.
    """
    has_trifluoromethyl_in_final = False
    has_trifluoromethyl_in_intermediate = False

    def dfs_traverse(node, depth=0):
        nonlocal has_trifluoromethyl_in_final, has_trifluoromethyl_in_intermediate

        if node["type"] == "mol":
            if "smiles" in node:
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol:
                    # Check for trifluoromethyl group
                    trifluoromethyl_pattern = Chem.MolFromSmarts("[#6]([#9])([#9])[#9]")
                    if mol.HasSubstructMatch(trifluoromethyl_pattern):
                        print(
                            f"Found trifluoromethyl group in molecule at depth {depth}: {node['smiles']}"
                        )
                        if depth == 0:  # Final product
                            has_trifluoromethyl_in_final = True
                        else:  # Intermediate
                            has_trifluoromethyl_in_intermediate = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Return True if trifluoromethyl is present in both final product and at least one intermediate
    strategy_present = (
        has_trifluoromethyl_in_final and has_trifluoromethyl_in_intermediate
    )
    print(f"Trifluoromethyl in final: {has_trifluoromethyl_in_final}")
    print(f"Trifluoromethyl in intermediate: {has_trifluoromethyl_in_intermediate}")
    print(f"Strategy detection result: {strategy_present}")
    return strategy_present
