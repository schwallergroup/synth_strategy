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
    This function detects if a nitro group is present throughout the synthesis without modification.
    """
    nitro_present_in_all_steps = True

    # SMARTS pattern for nitro group
    nitro_pattern = Chem.MolFromSmarts("[#7+](=[#8])[#8-]")

    def dfs_traverse(node):
        nonlocal nitro_present_in_all_steps

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            product_part = rsmi.split(">")[-1]
            product_mol = Chem.MolFromSmiles(product_part)

            if product_mol:
                if not product_mol.HasSubstructMatch(nitro_pattern):
                    nitro_present_in_all_steps = False
                    print(f"Nitro group not present in: {rsmi}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    print(f"Nitro group persistence strategy: {nitro_present_in_all_steps}")
    return nitro_present_in_all_steps
