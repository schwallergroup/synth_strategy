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
    Detects a linear synthesis strategy maintaining a fluorinated aromatic system throughout.
    """
    # Track the number of reactions and whether they're all linear
    reaction_count = 0
    is_linear = True
    maintains_fluorinated_aromatic = True

    def dfs_traverse(node):
        nonlocal reaction_count, is_linear, maintains_fluorinated_aromatic

        if node["type"] == "reaction":
            reaction_count += 1

            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check if more than one major fragment in product (non-linear)
            # Simplistic check: if product has more than one "." separator
            if product_smiles.count(".") > 0:
                is_linear = False

            # Check for fluorinated aromatic system
            product = Chem.MolFromSmiles(product_smiles)
            fluoro_aromatic_pattern = Chem.MolFromSmarts("[c][F]")

            if product is not None and not product.HasSubstructMatch(
                fluoro_aromatic_pattern
            ):
                maintains_fluorinated_aromatic = False

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if synthesis is linear with at least 3 steps and maintains fluorinated aromatic
    return is_linear and reaction_count >= 3 and maintains_fluorinated_aromatic
