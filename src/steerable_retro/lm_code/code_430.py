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
    Detects a linear strategy for elaborating a heterocyclic core through sequential reactions.
    """
    reaction_count = 0
    heterocycle_reactions = 0

    def dfs_traverse(node):
        nonlocal reaction_count, heterocycle_reactions

        if node["type"] == "reaction":
            reaction_count += 1

            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if reaction involves heterocycles
            heterocycle_in_reactants = False
            for reactant in reactants:
                mol = Chem.MolFromSmiles(reactant)
                if mol:
                    # Look for heterocyclic patterns (nitrogen-containing rings)
                    if mol.HasSubstructMatch(Chem.MolFromSmarts("[n]")) or mol.HasSubstructMatch(
                        Chem.MolFromSmarts("[N]1[C][C][C][C]1")
                    ):
                        heterocycle_in_reactants = True
                        break

            # Check if product also has heterocycle
            if heterocycle_in_reactants:
                prod_mol = Chem.MolFromSmiles(product)
                if prod_mol and (
                    prod_mol.HasSubstructMatch(Chem.MolFromSmarts("[n]"))
                    or prod_mol.HasSubstructMatch(Chem.MolFromSmarts("[N]1[C][C][C][C]1"))
                ):
                    heterocycle_reactions += 1
                    print(f"Heterocycle elaboration step detected (count: {heterocycle_reactions})")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Check if we have a linear synthesis (all reactions involve heterocycles)
    # and at least 4 reactions total
    return (
        reaction_count >= 4
        and heterocycle_reactions >= 3
        and heterocycle_reactions == reaction_count
    )
