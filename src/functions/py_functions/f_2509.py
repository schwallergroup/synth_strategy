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
    This function detects a strategy involving sequential SNAr reactions on a heterocyclic core.
    It looks for multiple reactions where a halogen on a heterocycle is replaced by O or N.
    """
    snar_reactions = 0

    def dfs_traverse(node):
        nonlocal snar_reactions

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for SNAr reaction pattern
            # Look for heterocycle with halogen in reactants and O/N connection in product
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
            product_mol = Chem.MolFromSmiles(product) if product else None

            if product_mol:
                # Pattern for heterocycle with halogen (Cl, F, Br, I)
                heterocycle_hal_pattern = Chem.MolFromSmarts(
                    "[c,n]1[c,n][c,n][c,n][c,n]1[Cl,F,Br,I]"
                )

                # Patterns for O/N connections in heterocycles
                o_connection_pattern = Chem.MolFromSmarts(
                    "[c,n]1[c,n][c,n][c,n][c,n]1[O]"
                )
                n_connection_pattern = Chem.MolFromSmarts(
                    "[c,n]1[c,n][c,n][c,n][c,n]1[N]"
                )

                # Check if any reactant has the halogen pattern
                has_hal_reactant = any(
                    mol and mol.HasSubstructMatch(heterocycle_hal_pattern)
                    for mol in reactant_mols
                )

                # Check if product has O or N connection
                has_o_product = product_mol.HasSubstructMatch(o_connection_pattern)
                has_n_product = product_mol.HasSubstructMatch(n_connection_pattern)

                if has_hal_reactant and (has_o_product or has_n_product):
                    snar_reactions += 1
                    print(f"SNAr reaction detected: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if at least 2 SNAr reactions are detected
    return snar_reactions >= 2
