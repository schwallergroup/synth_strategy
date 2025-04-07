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
    Detects if the synthesis route employs a linear chain extension strategy
    with a 3-carbon propyl linker
    """
    chain_extension_found = False

    def dfs_traverse(node):
        nonlocal chain_extension_found

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Look for dibromoalkane pattern
            dibromo_pattern = Chem.MolFromSmarts("[Br][C][C][C][Br]")

            # Look for aryl ether with propyl chain pattern in product
            aryl_ether_pattern = Chem.MolFromSmarts("[c][O][C][C][C][Br]")

            dibromo_found = False
            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if reactant_mol and reactant_mol.HasSubstructMatch(dibromo_pattern):
                    dibromo_found = True
                    break

            product_mol = Chem.MolFromSmiles(product)
            if (
                dibromo_found
                and product_mol
                and product_mol.HasSubstructMatch(aryl_ether_pattern)
            ):
                print("Found linear chain extension with 3-carbon propyl linker")
                chain_extension_found = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return chain_extension_found
