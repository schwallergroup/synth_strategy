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
    This function detects if the synthesis involves an intramolecular cyclization
    (number of fragments decreases in a reaction).
    """
    has_intramolecular_cyclization = False

    def dfs_traverse(node, depth=0):
        nonlocal has_intramolecular_cyclization

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Count number of fragments
            reactant_fragments = len(reactants_smiles)
            product_fragments = len(product_smiles.split("."))

            # If product has fewer fragments than reactants, it's an intramolecular reaction
            if product_fragments < reactant_fragments and product_fragments == 1:
                print(f"Intramolecular cyclization detected at depth {depth}")
                has_intramolecular_cyclization = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    return has_intramolecular_cyclization
