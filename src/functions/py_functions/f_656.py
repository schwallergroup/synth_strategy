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
    This function detects multiple SNAr reactions in the synthetic route.
    """
    snar_reactions = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            # Extract reactants and products
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0]
            product = rsmi.split(">")[-1]

            # Check for SNAr pattern: ArX + nucleophile → Ar-nucleophile
            # Look for halogen displacement by O, N, or S nucleophiles
            reactants_mol = Chem.MolFromSmiles(reactants)
            product_mol = Chem.MolFromSmiles(product)

            if reactants_mol and product_mol:
                # Check for aryl halide in reactants
                aryl_halide_pattern = Chem.MolFromSmarts("[c]-[F,Cl,Br,I]")

                # Check for aryl ether, amine, or thioether in product
                aryl_o_pattern = Chem.MolFromSmarts("[c]-[O]")
                aryl_n_pattern = Chem.MolFromSmarts("[c]-[N]")
                aryl_s_pattern = Chem.MolFromSmarts("[c]-[S]")

                if reactants_mol.HasSubstructMatch(aryl_halide_pattern) and (
                    product_mol.HasSubstructMatch(aryl_o_pattern)
                    or product_mol.HasSubstructMatch(aryl_n_pattern)
                    or product_mol.HasSubstructMatch(aryl_s_pattern)
                ):
                    snar_reactions.append(depth)
                    print(f"SNAr reaction found at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if at least 2 SNAr reactions were found
    if len(snar_reactions) >= 2:
        print(f"Sequential SNAr strategy detected with {len(snar_reactions)} reactions")
        return True
    return False
