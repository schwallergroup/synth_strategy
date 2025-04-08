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
    This function detects a synthetic strategy involving ester hydrolysis to form
    a carboxylic acid.
    """
    # Track if we found ester hydrolysis
    found_ester_hydrolysis = False

    def dfs_traverse(node):
        nonlocal found_ester_hydrolysis

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for ester to carboxylic acid conversion
            ester_pattern = Chem.MolFromSmarts("[C$(C=O)][O][C]")
            carboxylic_acid_pattern = Chem.MolFromSmarts("[C$(C=O)][OH]")

            # Check reactants for ester
            has_ester = False
            for reactant in reactants:
                mol = Chem.MolFromSmiles(reactant)
                if mol and mol.HasSubstructMatch(ester_pattern):
                    has_ester = True
                    break

            # Check product for carboxylic acid
            product_mol = Chem.MolFromSmiles(product)
            has_carboxylic_acid = product_mol and product_mol.HasSubstructMatch(
                carboxylic_acid_pattern
            )

            if has_ester and has_carboxylic_acid:
                found_ester_hydrolysis = True
                print("Found ester hydrolysis strategy")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return found_ester_hydrolysis
