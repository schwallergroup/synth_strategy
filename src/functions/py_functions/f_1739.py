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
    This function detects a synthetic strategy involving the construction of
    a nitrogen-containing heterocyclic system.
    """
    n_heterocycle_formed = False

    def dfs_traverse(node, depth=0):
        nonlocal n_heterocycle_formed

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for nitrogen heterocycle formation
            product_mol = Chem.MolFromSmiles(product)
            reactants_mols = [
                Chem.MolFromSmiles(r) for r in reactants if Chem.MolFromSmiles(r)
            ]

            if product_mol and reactants_mols:
                # Check for nitrogen in aromatic ring in product
                n_heterocycle_pattern = Chem.MolFromSmarts("[n;R]")
                if product_mol.HasSubstructMatch(n_heterocycle_pattern):
                    # Check if any reactant doesn't have this pattern
                    if any(
                        not r.HasSubstructMatch(n_heterocycle_pattern)
                        for r in reactants_mols
                    ):
                        n_heterocycle_formed = True
                        print(
                            f"Nitrogen heterocycle formation detected at depth {depth}"
                        )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    if n_heterocycle_formed:
        print("Nitrogen heterocycle construction strategy detected")

    return n_heterocycle_formed
