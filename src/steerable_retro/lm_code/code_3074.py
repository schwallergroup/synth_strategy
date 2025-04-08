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
    Detects the use of protection/deprotection strategies,
    particularly Boc protection of nitrogen.
    """
    protection_found = False

    def dfs_traverse(node):
        nonlocal protection_found

        if node["type"] == "reaction":
            rsmi = node["metadata"].get("rsmi", "")
            if not rsmi:
                return

            product_smiles = rsmi.split(">")[-1]

            try:
                product_mol = Chem.MolFromSmiles(product_smiles)
                if product_mol is None:
                    return

                # Check for Boc protected amine
                boc_pattern = Chem.MolFromSmarts("[#6]-[#8]-[#6](=[#8])-[#7]")

                if product_mol.HasSubstructMatch(boc_pattern):
                    protection_found = True
                    print("Found Boc protection")

            except:
                pass  # Skip if there's an error processing the molecules

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Protection/deprotection strategy detected: {protection_found}")
    return protection_found
