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
    This function detects a synthetic strategy involving benzyl protection/deprotection
    of an alcohol.
    """
    # Track if we found benzyl deprotection
    found_benzyl_deprotection = False

    def dfs_traverse(node):
        nonlocal found_benzyl_deprotection

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for benzyl ether to alcohol conversion (deprotection)
            benzyl_ether_pattern = Chem.MolFromSmarts("[O][C][c]1[c][c][c][c][c]1")
            alcohol_pattern = Chem.MolFromSmarts("[OH]")

            # Check reactants for benzyl ether
            has_benzyl_ether = False
            for reactant in reactants:
                mol = Chem.MolFromSmiles(reactant)
                if mol and mol.HasSubstructMatch(benzyl_ether_pattern):
                    has_benzyl_ether = True
                    break

            # Check product for alcohol
            product_mol = Chem.MolFromSmiles(product)
            has_alcohol = product_mol and product_mol.HasSubstructMatch(alcohol_pattern)

            if has_benzyl_ether and has_alcohol:
                found_benzyl_deprotection = True
                print("Found benzyl deprotection strategy")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return found_benzyl_deprotection
