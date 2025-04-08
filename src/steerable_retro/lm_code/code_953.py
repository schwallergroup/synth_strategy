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
    Detects if the synthesis route includes a halogen exchange step
    (bromide to iodide conversion)
    """
    halogen_exchange_found = False

    def dfs_traverse(node):
        nonlocal halogen_exchange_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Look for C-Br pattern in reactants
            c_br_pattern = Chem.MolFromSmarts("[C][Br]")
            # Look for C-I pattern in product
            c_i_pattern = Chem.MolFromSmarts("[C][I]")

            reactant_mol = Chem.MolFromSmiles(reactants[0])
            product_mol = Chem.MolFromSmiles(product)

            if reactant_mol and product_mol:
                if reactant_mol.HasSubstructMatch(c_br_pattern) and product_mol.HasSubstructMatch(
                    c_i_pattern
                ):
                    print("Found halogen exchange (Br to I)")
                    halogen_exchange_found = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return halogen_exchange_found
