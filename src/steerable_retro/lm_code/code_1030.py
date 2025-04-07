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
    Detects if the synthesis includes an epoxide opening reaction.
    """
    found_epoxide_opening = False

    def dfs_traverse(node):
        nonlocal found_epoxide_opening

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Convert to RDKit molecules
            reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles]
            product = Chem.MolFromSmiles(product_smiles)

            # Check for epoxide opening
            epoxide_pattern = Chem.MolFromSmarts("[C]1[O][C]1")
            opened_pattern = Chem.MolFromSmarts("[C][C][O]")

            if any(
                r.HasSubstructMatch(epoxide_pattern) for r in reactants
            ) and product.HasSubstructMatch(opened_pattern):
                found_epoxide_opening = True
                print("Found epoxide opening")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return found_epoxide_opening
