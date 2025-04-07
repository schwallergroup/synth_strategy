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
    Detects synthesis strategy involving difluoromethoxy ether formation
    from a phenol group.
    """
    difluoromethoxy_formation = False

    def dfs_traverse(node):
        nonlocal difluoromethoxy_formation

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if any reactant has a phenol group
            phenol_present = any(
                Chem.MolFromSmiles(r)
                and Chem.MolFromSmiles(r).HasSubstructMatch(
                    Chem.MolFromSmarts("[c][OH]")
                )
                for r in reactants
                if Chem.MolFromSmiles(r)
            )

            # Check if product has a difluoromethoxy group
            product_mol = Chem.MolFromSmiles(product)
            if product_mol and product_mol.HasSubstructMatch(
                Chem.MolFromSmarts("[c][O][CH]([F])[F]")
            ):
                if phenol_present:
                    difluoromethoxy_formation = True
                    print("Detected difluoromethoxy ether formation from phenol")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(
        f"Difluoromethoxy ether formation strategy detected: {difluoromethoxy_formation}"
    )
    return difluoromethoxy_formation
