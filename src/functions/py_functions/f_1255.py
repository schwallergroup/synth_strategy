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
    This function detects a synthetic strategy involving piperazine ring formation.
    """
    has_piperazine_formation = False

    def dfs_traverse(node):
        nonlocal has_piperazine_formation

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"].get("rsmi", "")
            if rsmi:
                reactants_part = rsmi.split(">")[0]
                product_part = rsmi.split(">")[-1]

                reactants = [
                    Chem.MolFromSmiles(r) for r in reactants_part.split(".") if r
                ]
                product = Chem.MolFromSmiles(product_part) if product_part else None

                if product and all(r for r in reactants):
                    # Check for piperazine in product
                    piperazine_pattern = Chem.MolFromSmarts(
                        "[#7]1[#6][#6][#7][#6][#6]1"
                    )

                    # Check for bis-chloroethylamine pattern in reactants
                    bis_chloro_pattern = Chem.MolFromSmarts("Cl[#6][#6][#7][#6][#6]Cl")

                    has_piperazine = product.HasSubstructMatch(piperazine_pattern)
                    has_bis_chloro = any(
                        r.HasSubstructMatch(bis_chloro_pattern) for r in reactants
                    )

                    # Check for aniline in reactants
                    aniline_pattern = Chem.MolFromSmarts("[c][NH2]")
                    has_aniline = any(
                        r.HasSubstructMatch(aniline_pattern) for r in reactants
                    )

                    if has_piperazine and (has_bis_chloro or has_aniline):
                        print("Detected piperazine ring formation")
                        has_piperazine_formation = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Piperazine formation strategy detected: {has_piperazine_formation}")
    return has_piperazine_formation
