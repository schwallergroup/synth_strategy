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
    Detects if the synthesis route preserves a trifluoromethyl group throughout the synthesis.
    """
    final_product_has_cf3 = False
    starting_material_has_cf3 = False

    def dfs_traverse(node):
        nonlocal final_product_has_cf3, starting_material_has_cf3

        if node["type"] == "mol":
            cf3_pattern = Chem.MolFromSmarts("[C]([F])([F])[F]")
            mol = Chem.MolFromSmiles(node["smiles"])

            if mol and mol.HasSubstructMatch(cf3_pattern):
                # Check if this is the final product
                if "depth" in node and node["depth"] == 0:
                    final_product_has_cf3 = True
                    print("Final product contains trifluoromethyl group")

                # Check if this is a starting material
                elif node.get("in_stock", False):
                    starting_material_has_cf3 = True
                    print("Starting material contains trifluoromethyl group")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if trifluoromethyl group is found in both starting materials and final product
    return final_product_has_cf3 and starting_material_has_cf3
