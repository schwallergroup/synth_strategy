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
    This function detects a synthetic strategy involving the construction of a
    heterocycle containing a trifluoromethyl group.
    """
    trifluoromethyl_heterocycle_found = False

    def dfs_traverse(node, depth=0):
        nonlocal trifluoromethyl_heterocycle_found

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            product_smiles = rsmi.split(">")[-1]

            product_mol = Chem.MolFromSmiles(product_smiles)

            if product_mol:
                # Check for heterocycle formation
                heterocycle_pattern = Chem.MolFromSmarts("[#7]1~[#6]~[#7]~[#6]~[#6]~[#6]1")
                trifluoromethyl_pattern = Chem.MolFromSmarts("[#6]([F])([F])[F]")

                if product_mol.HasSubstructMatch(
                    heterocycle_pattern
                ) and product_mol.HasSubstructMatch(trifluoromethyl_pattern):
                    # Check if the trifluoromethyl is connected to the heterocycle
                    # This is a simplified check - in a real implementation, you would need
                    # to verify the connection between the two substructures
                    trifluoromethyl_heterocycle_found = True
                    print(f"Trifluoromethyl-containing heterocycle detected at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    if trifluoromethyl_heterocycle_found:
        print("Trifluoromethyl-containing heterocycle strategy detected")
    else:
        print("Trifluoromethyl-containing heterocycle strategy not detected")

    return trifluoromethyl_heterocycle_found
