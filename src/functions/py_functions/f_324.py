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
    This function detects nitro reduction (NO2 -> NH2) in the synthetic route.
    """
    nitro_reduction_found = False

    def dfs_traverse(node):
        nonlocal nitro_reduction_found

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0]
            product = rsmi.split(">")[-1]

            # Check if reactants contain nitro group
            reactants_mol = Chem.MolFromSmiles(reactants)
            if (
                reactants_mol
                and reactants_mol.HasSubstructMatch(Chem.MolFromSmarts("[N+](=O)[O-]"))
                or reactants_mol.HasSubstructMatch(Chem.MolFromSmarts("[N](=O)(=O)"))
            ):

                # Check if product contains amine group where nitro was
                product_mol = Chem.MolFromSmiles(product)
                if product_mol and product_mol.HasSubstructMatch(
                    Chem.MolFromSmarts("[NH2]")
                ):
                    print("Nitro reduction detected")
                    nitro_reduction_found = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return nitro_reduction_found
