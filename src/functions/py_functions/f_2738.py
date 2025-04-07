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
    This function detects a synthetic strategy involving conversion of an alcohol
    to a chloride.
    """
    # Track if we found alcohol to chloride conversion
    found_alcohol_to_chloride = False

    def dfs_traverse(node):
        nonlocal found_alcohol_to_chloride

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for alcohol to chloride conversion
            alcohol_pattern = Chem.MolFromSmarts("[OH]")
            chloride_pattern = Chem.MolFromSmarts("[Cl]")

            # Check reactants for alcohol
            has_alcohol = False
            for reactant in reactants:
                mol = Chem.MolFromSmiles(reactant)
                if mol and mol.HasSubstructMatch(alcohol_pattern):
                    has_alcohol = True
                    break

            # Check product for chloride
            product_mol = Chem.MolFromSmiles(product)
            has_chloride = product_mol and product_mol.HasSubstructMatch(
                chloride_pattern
            )

            if has_alcohol and has_chloride:
                found_alcohol_to_chloride = True
                print("Found alcohol to chloride conversion strategy")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return found_alcohol_to_chloride
