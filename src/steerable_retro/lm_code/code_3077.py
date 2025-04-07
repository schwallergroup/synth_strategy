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
    Detects a synthetic route that includes reduction of a nitro group to an amine.
    """
    has_nitro_reduction = False

    def dfs_traverse(node):
        nonlocal has_nitro_reduction

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_str = rsmi.split(">")[0]
                product_str = rsmi.split(">")[-1]

                # Check for nitro group in reactants
                nitro_pattern = Chem.MolFromSmarts("[#6]-[N+](=[O])[O-]")
                aniline_pattern = Chem.MolFromSmarts("c-[NH2]")

                reactants = [Chem.MolFromSmiles(r) for r in reactants_str.split(".") if r]
                product = Chem.MolFromSmiles(product_str) if product_str else None

                if product and any(
                    r and r.HasSubstructMatch(nitro_pattern) for r in reactants if r
                ):
                    if product.HasSubstructMatch(aniline_pattern):
                        has_nitro_reduction = True
                        print("Found nitro reduction to aniline")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Nitro reduction strategy detection result: {has_nitro_reduction}")
    return has_nitro_reduction
