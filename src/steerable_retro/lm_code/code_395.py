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
    Detects a strategy involving late-stage amination (Cl → NH₂ transformation).
    """
    # Track if we found late-stage amination
    amination_found = False

    def dfs_traverse(node, depth=0):
        nonlocal amination_found

        if node["type"] == "reaction" and depth <= 1:  # Consider only late-stage reactions
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            try:
                # Check for Cl → NH₂ transformation
                if any("Cl" in r for r in reactants_smiles) and "NH2" in product_smiles:
                    amination_found = True
                    print(f"Late-stage amination detected at depth {depth}")

            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    if amination_found:
        print("Late-stage amination strategy detected")

    return amination_found
