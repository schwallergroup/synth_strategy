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
    Detects if the synthesis preserves stereochemistry from starting materials.
    """
    has_stereochemistry = False
    preserves_stereochemistry = True

    def dfs_traverse(node, depth=0):
        nonlocal has_stereochemistry, preserves_stereochemistry

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for stereochemistry in SMILES
            stereo_pattern = re.compile(r"[@]")

            # Check product for stereochemistry
            if stereo_pattern.search(product_smiles):
                has_stereochemistry = True
                print(f"Found stereochemistry in product at depth {depth}")

            # Check if stereochemistry is preserved from reactants
            reactants_have_stereo = any(
                stereo_pattern.search(r) for r in reactants_smiles
            )
            product_has_stereo = stereo_pattern.search(product_smiles)

            if reactants_have_stereo and not product_has_stereo:
                preserves_stereochemistry = False
                print(f"Stereochemistry lost at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Return True if stereochemistry is present and preserved
    return has_stereochemistry and preserves_stereochemistry
