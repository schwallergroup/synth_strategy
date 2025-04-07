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
    This function detects a strategy involving modification of nitrogen-containing heterocycles.
    """
    heterocycle_modifications = 0

    def dfs_traverse(node, depth=0):
        nonlocal heterocycle_modifications

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]

            try:
                # Check for heterocycle modification
                product_mol = Chem.MolFromSmiles(product_part)

                # Look for nitrogen heterocycles in product
                n_heterocycle_pattern = Chem.MolFromSmarts(
                    "[n]1[c,n][c,n][c,n][c,n][c,n]1"
                )  # 6-membered N-heterocycle

                if product_mol and product_mol.HasSubstructMatch(n_heterocycle_pattern):
                    # Check if any reactant also has this pattern
                    reactants = reactants_part.split(".")
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol and reactant_mol.HasSubstructMatch(
                            n_heterocycle_pattern
                        ):
                            # The heterocycle exists in both reactant and product, suggesting modification
                            heterocycle_modifications += 1
                            print(f"Found heterocycle modification at depth {depth}")
                            break
            except:
                print(f"Error processing reaction at depth {depth}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Return True if at least one heterocycle modification was found
    return heterocycle_modifications > 0
