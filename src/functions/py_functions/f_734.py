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
    Detects synthesis strategy involving reduction of aldehyde to alcohol.
    """
    aldehyde_reduction_found = False

    def dfs_traverse(node):
        nonlocal aldehyde_reduction_found

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for aldehyde to alcohol reduction
            aldehyde_pattern = Chem.MolFromSmarts("[C:1]=[O:2]")
            alcohol_pattern = Chem.MolFromSmarts("[C:1]-[O:2]")

            product_mol = Chem.MolFromSmiles(product)
            if product_mol and product_mol.HasSubstructMatch(alcohol_pattern):
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(
                        aldehyde_pattern
                    ):
                        print("Detected aldehyde reduction to alcohol")
                        aldehyde_reduction_found = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return aldehyde_reduction_found
