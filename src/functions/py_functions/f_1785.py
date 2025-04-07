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
    This function detects late-stage carbamate formation.
    """
    carbamate_formation_found = False

    def dfs_traverse(node, depth=0):
        nonlocal carbamate_formation_found

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]

            # Check if product has carbamate
            product_mol = Chem.MolFromSmiles(product_part)

            if product_mol:
                carbamate_pattern = Chem.MolFromSmarts("[N][C](=[O])[O]")
                if product_mol.HasSubstructMatch(carbamate_pattern):
                    # Check if this is a late-stage reaction (depth < 2)
                    if depth < 2:
                        print(f"Found late-stage carbamate formation at depth {depth}")
                        carbamate_formation_found = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return carbamate_formation_found
