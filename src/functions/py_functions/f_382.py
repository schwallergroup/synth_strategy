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
    Detects early-stage thioether (C-S-C) formation.
    """
    found_thioether_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal found_thioether_formation

        if (
            node["type"] == "reaction" and depth >= 3
        ):  # Focus on early-stage reactions (high depth)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for thiol pattern in reactants
                thiol_pattern = Chem.MolFromSmarts("[#16;H1]")

                # Check for thioether pattern in product
                thioether_pattern = Chem.MolFromSmarts("[#6]-[#16]-[#6]")

                has_thiol = any(
                    Chem.MolFromSmiles(r)
                    and Chem.MolFromSmiles(r).HasSubstructMatch(thiol_pattern)
                    for r in reactants
                    if r
                )

                product_mol = Chem.MolFromSmiles(product)
                has_thioether_in_product = (
                    product_mol and product_mol.HasSubstructMatch(thioether_pattern)
                )

                if has_thiol and has_thioether_in_product:
                    found_thioether_formation = True
                    print(f"Found early-stage thioether formation at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    return found_thioether_formation
