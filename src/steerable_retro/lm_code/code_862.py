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
    Detects a synthetic route that includes use of protecting groups (specifically TIPS).
    """
    has_protecting_group = False

    def dfs_traverse(node):
        nonlocal has_protecting_group

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for silicon-containing protecting groups
            if "Si" in "".join(reactants_smiles) or "Si" in product_smiles:
                has_protecting_group = True
                print("Found protecting group (Si-containing)")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(f"Protecting group strategy: {has_protecting_group}")
    return has_protecting_group
