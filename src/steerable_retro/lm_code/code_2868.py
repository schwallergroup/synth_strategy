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
    Detects if the synthetic route includes a nitro reduction step.
    Looks for transformation of nitro group to amine.
    """
    found_nitro_reduction = False

    def is_nitro_reduction(reactants, product):
        # Check if reactant has nitro group
        nitro_pattern = Chem.MolFromSmarts("[#6]-[N+](=[O])[O-]")
        has_nitro = False
        for r in reactants:
            r_mol = Chem.MolFromSmiles(r)
            if r_mol and r_mol.HasSubstructMatch(nitro_pattern):
                has_nitro = True
                break

        if not has_nitro:
            return False

        # Check if product has amine group where nitro was
        amine_pattern = Chem.MolFromSmarts("[#6]-[NH2]")
        product_mol = Chem.MolFromSmiles(product)
        if product_mol and product_mol.HasSubstructMatch(amine_pattern):
            return True

        return False

    def dfs_traverse(node):
        nonlocal found_nitro_reduction

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            if is_nitro_reduction(reactants, product):
                found_nitro_reduction = True
                print("Found nitro reduction")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    print(f"Nitro reduction strategy: {found_nitro_reduction}")
    return found_nitro_reduction
