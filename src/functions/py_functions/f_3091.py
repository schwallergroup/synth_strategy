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
    This function detects a strategy involving amide bond formation.
    """
    amide_formation_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal amide_formation_detected

        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                product_part = rsmi.split(">")[-1]

                # Check for amide bond in product
                amide_pattern = Chem.MolFromSmarts("[#7]-C(=O)-[#6]")
                product_mol = Chem.MolFromSmiles(product_part)

                if product_mol and product_mol.HasSubstructMatch(amide_pattern):
                    # Check if this is a new amide formation (simplified approach)
                    print(f"Found amide bond formation at depth {depth}")
                    amide_formation_detected = True

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    return amide_formation_detected
