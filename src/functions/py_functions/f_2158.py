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
    This function detects if a fluorinated aromatic is maintained throughout the synthesis.
    """
    fluoro_aromatic_pattern = Chem.MolFromSmarts("[c][F]")
    reactions_with_fluoro = 0
    total_reactions = 0

    def dfs_traverse(node):
        nonlocal reactions_with_fluoro, total_reactions

        if node["type"] == "reaction":
            total_reactions += 1
            try:
                rsmi = node["metadata"]["rsmi"]
                product_smiles = rsmi.split(">")[-1]

                product_mol = Chem.MolFromSmiles(product_smiles)
                if product_mol and product_mol.HasSubstructMatch(
                    fluoro_aromatic_pattern
                ):
                    reactions_with_fluoro += 1
                    print(
                        f"Fluorinated aromatic detected in reaction product: {product_smiles}"
                    )
            except Exception as e:
                print(f"Error in fluorinated aromatic detection: {e}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    # Consider fluorinated aromatic maintained if present in at least 75% of reactions
    return total_reactions > 0 and reactions_with_fluoro / total_reactions >= 0.75
