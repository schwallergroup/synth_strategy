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
    This function detects if a benzylic alcohol motif is maintained
    throughout most of the synthesis.
    """
    reactions_with_benzylic_alcohol = 0
    total_reactions = 0

    def dfs_traverse(node):
        nonlocal reactions_with_benzylic_alcohol, total_reactions

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                total_reactions += 1
                rsmi = node["metadata"]["rsmi"]
                product = rsmi.split(">")[-1]

                # Check for benzylic alcohol pattern
                benzylic_alcohol_pattern = Chem.MolFromSmarts("[c]-[CH](-[OH])-[#6]")
                product_mol = Chem.MolFromSmiles(product)

                if product_mol and product_mol.HasSubstructMatch(
                    benzylic_alcohol_pattern
                ):
                    reactions_with_benzylic_alcohol += 1
                    print("Benzylic alcohol detected in reaction product")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from root
    dfs_traverse(route)

    # If benzylic alcohol is present in at least 60% of reactions
    has_benzylic_alcohol_throughout = (
        total_reactions > 0
        and (reactions_with_benzylic_alcohol / total_reactions) >= 0.6
    )

    if has_benzylic_alcohol_throughout:
        print(
            f"Benzylic alcohol maintained throughout synthesis: {reactions_with_benzylic_alcohol}/{total_reactions} reactions"
        )

    return has_benzylic_alcohol_throughout
