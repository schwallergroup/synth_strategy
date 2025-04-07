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
    This function detects a synthetic strategy involving nitro reduction to form an amine.
    """
    has_nitro_reduction = False

    def dfs_traverse(node):
        nonlocal has_nitro_reduction

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"].get("rsmi", "")
            if rsmi:
                reactants_part = rsmi.split(">")[0]
                product_part = rsmi.split(">")[-1]

                reactants = [Chem.MolFromSmiles(r) for r in reactants_part.split(".") if r]
                product = Chem.MolFromSmiles(product_part) if product_part else None

                if product and all(r for r in reactants):
                    # Check for nitro group in reactant
                    nitro_pattern = Chem.MolFromSmarts("[c][N+](=[O])[O-]")
                    # Check for amine in product
                    amine_pattern = Chem.MolFromSmarts("[c][NH2]")

                    has_nitro = any(r.HasSubstructMatch(nitro_pattern) for r in reactants)
                    has_amine = product.HasSubstructMatch(amine_pattern)

                    if has_nitro and has_amine:
                        print("Detected nitro reduction to amine")
                        has_nitro_reduction = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Nitro reduction strategy detected: {has_nitro_reduction}")
    return has_nitro_reduction
