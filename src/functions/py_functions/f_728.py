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
    This function detects a linear synthesis strategy using a propyl linker
    to connect fragments.
    """
    propyl_linker = False
    linear_synthesis = True  # Assume linear until proven otherwise
    reaction_count = 0

    def dfs_traverse(node):
        nonlocal propyl_linker, linear_synthesis, reaction_count

        if node["type"] == "reaction":
            reaction_count += 1

            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                product = rsmi.split(">")[-1]

                # Check for propyl linker between nitrogens
                product_mol = Chem.MolFromSmiles(product)
                if product_mol:
                    propyl_linker_pattern = Chem.MolFromSmarts("[N][C][C][C][N]")
                    if product_mol.HasSubstructMatch(propyl_linker_pattern):
                        propyl_linker = True
                        print("Detected propyl linker between nitrogens")

            # Check if this is a convergent step (more than one non-trivial reactant)
            if "children" in node:
                non_trivial_children = 0
                for child in node["children"]:
                    if child["type"] == "mol" and not child.get("in_stock", False):
                        non_trivial_children += 1

                if non_trivial_children > 1:
                    linear_synthesis = False
                    print("Detected convergent synthesis step")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return propyl_linker and linear_synthesis and reaction_count >= 3
