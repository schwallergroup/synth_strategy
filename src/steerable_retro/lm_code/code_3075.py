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
    Detects a synthetic route with multiple SNAr reactions (at least 2).
    """
    snar_count = 0

    def dfs_traverse(node):
        nonlocal snar_count

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_str = rsmi.split(">")[0]
                product_str = rsmi.split(">")[-1]

                reactants = [Chem.MolFromSmiles(r) for r in reactants_str.split(".") if r]
                product = Chem.MolFromSmiles(product_str) if product_str else None

                if product and reactants:
                    # Check for SNAr pattern: aryl halide + amine
                    halide_pattern = Chem.MolFromSmarts("[c]-[F,Cl,Br,I]")
                    amine_pattern = Chem.MolFromSmarts("[N;H0,H1,H2]")

                    has_halide = any(
                        r and r.HasSubstructMatch(halide_pattern) for r in reactants if r
                    )
                    has_amine = any(
                        r and r.HasSubstructMatch(amine_pattern) for r in reactants if r
                    )

                    # Check if product has new C-N bond where halide was
                    if has_halide and has_amine:
                        snar_count += 1
                        print(f"Found potential SNAr reaction, count: {snar_count}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if at least 2 SNAr reactions are found
    result = snar_count >= 2
    print(f"Multiple SNAr strategy detection result: {result}")
    return result
