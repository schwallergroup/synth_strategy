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
    Detects a strategy involving reduction of an aromatic nitro group to an amine.
    """
    nitro_reduction_found = False

    def dfs_traverse(node):
        nonlocal nitro_reduction_found

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for nitro group in reactant
            reactant_has_nitro = False
            for reactant in reactants:
                mol = Chem.MolFromSmiles(reactant)
                if mol and mol.HasSubstructMatch(Chem.MolFromSmarts("[c][N+](=[O])[O-]")):
                    reactant_has_nitro = True
                    break

            # Check for amine in product where nitro was
            if reactant_has_nitro:
                prod_mol = Chem.MolFromSmiles(product)
                if prod_mol and prod_mol.HasSubstructMatch(Chem.MolFromSmarts("[c][NH2]")):
                    nitro_reduction_found = True
                    print("Nitro reduction detected")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return nitro_reduction_found
