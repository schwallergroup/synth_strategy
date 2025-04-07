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
    Detects if the synthesis includes N-alkylation of a heterocycle.
    """
    found_n_alkylation = False

    def dfs_traverse(node):
        nonlocal found_n_alkylation

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Convert to RDKit molecules
            reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles]
            product = Chem.MolFromSmiles(product_smiles)

            # Check for N-alkylation of heterocycle
            heterocycle_nh = Chem.MolFromSmarts("[nH]")
            alkyl_halide = Chem.MolFromSmarts("[C][Br,Cl,I,F]")
            n_alkylated = Chem.MolFromSmarts("[n][C]")

            if (
                any(r.HasSubstructMatch(heterocycle_nh) for r in reactants)
                and any(r.HasSubstructMatch(alkyl_halide) for r in reactants)
                and product.HasSubstructMatch(n_alkylated)
            ):
                found_n_alkylation = True
                print("Found N-alkylation of heterocycle")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return found_n_alkylation
