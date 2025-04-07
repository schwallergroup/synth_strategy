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
    This function detects if the synthesis route includes a Boc protection step.
    """
    boc_protection_found = False

    def dfs_traverse(node):
        nonlocal boc_protection_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if any reactant has an amine that becomes Boc-protected in the product
            amine_pattern = Chem.MolFromSmarts("[NH2]")
            boc_pattern = Chem.MolFromSmarts("[NH]C(=O)OC(C)(C)C")

            try:
                product_mol = Chem.MolFromSmiles(product)
                if product_mol and product_mol.HasSubstructMatch(boc_pattern):
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol and reactant_mol.HasSubstructMatch(amine_pattern):
                            print("Boc protection detected")
                            boc_protection_found = True
                            break
            except:
                pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return boc_protection_found
