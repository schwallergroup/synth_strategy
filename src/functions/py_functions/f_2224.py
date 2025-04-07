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
    This function detects a pyrazole ring formation strategy in the synthesis route.
    """
    # Track if we found the strategy
    found_strategy = False

    def dfs_traverse(node, depth=0):
        nonlocal found_strategy

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]

            # Check if product contains a pyrazole ring
            product_mol = Chem.MolFromSmiles(product_part)
            if product_mol:
                pyrazole_pattern = Chem.MolFromSmarts("[n]1[n][c][c][c]1")
                if product_mol.HasSubstructMatch(pyrazole_pattern):
                    # Check if reactants don't have the pyrazole pattern
                    reactants_have_pyrazole = False
                    for reactant in reactants_part.split("."):
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol and reactant_mol.HasSubstructMatch(
                            pyrazole_pattern
                        ):
                            reactants_have_pyrazole = True
                            break

                    if not reactants_have_pyrazole:
                        print(f"Found pyrazole formation at depth {depth}")
                        found_strategy = True

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return found_strategy
