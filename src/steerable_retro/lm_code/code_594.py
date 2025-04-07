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
    This function detects if the synthesis includes activation of a hydroxyl group
    to a better leaving group (e.g., mesylate).
    """
    has_leaving_group_activation = False

    def dfs_traverse(node):
        nonlocal has_leaving_group_activation

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for alcohol in reactants
            alcohol_pattern = Chem.MolFromSmarts("[OH][#6]")

            # Check for mesylate in product
            mesylate_pattern = Chem.MolFromSmarts("[#6]OS(=O)(=O)C")

            has_alcohol = False
            for reactant in reactants:
                mol = Chem.MolFromSmiles(reactant)
                if mol and mol.HasSubstructMatch(alcohol_pattern):
                    has_alcohol = True

            product_mol = Chem.MolFromSmiles(product)
            has_mesylate = product_mol and product_mol.HasSubstructMatch(mesylate_pattern)

            if has_alcohol and has_mesylate:
                has_leaving_group_activation = True
                print("Detected leaving group activation strategy")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return has_leaving_group_activation
