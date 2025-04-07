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
    This function detects a strategy based on maintaining an indole-benzamide core
    throughout the synthesis.
    """
    has_indole_benzamide = False
    indole_benzamide_count = 0
    total_reactions = 0

    def dfs_traverse(node, depth=0):
        nonlocal has_indole_benzamide, indole_benzamide_count, total_reactions

        if node["type"] == "reaction":
            total_reactions += 1
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                product = rsmi.split(">")[-1]

                # Check for indole-benzamide core
                prod_mol = Chem.MolFromSmiles(product)
                if prod_mol:
                    indole_pattern = Chem.MolFromSmarts("c1cnc2ccccc12")
                    benzamide_pattern = Chem.MolFromSmarts("c1ccccc1C(=O)N")

                    if prod_mol.HasSubstructMatch(
                        indole_pattern
                    ) and prod_mol.HasSubstructMatch(benzamide_pattern):
                        indole_benzamide_count += 1
                        print(f"Found indole-benzamide core at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Consider it a core strategy if present in most reactions
    has_indole_benzamide = indole_benzamide_count >= total_reactions * 0.7

    if has_indole_benzamide:
        print(
            f"Detected indole-benzamide core strategy (present in {indole_benzamide_count}/{total_reactions} reactions)"
        )

    return has_indole_benzamide
