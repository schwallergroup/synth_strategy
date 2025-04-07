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
    Detects if the synthesis route involves a late-stage thioether formation.
    Late stage is defined as occurring in the first half of the synthesis (lower depth numbers).
    """
    thioether_formed = False
    max_depth = 0
    thioether_depth = None

    def dfs_traverse(node, depth=0):
        nonlocal thioether_formed, max_depth, thioether_depth

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this reaction forms a thioether bond (C-S-C)
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
            product_mol = Chem.MolFromSmiles(product) if product else None

            if product_mol:
                # Check for thioether pattern in product but not in reactants
                thioether_pattern = Chem.MolFromSmarts("[#6]-[#16]-[#6]")
                if product_mol.HasSubstructMatch(thioether_pattern):
                    # Check if thioether was not present in reactants
                    reactants_have_thioether = any(
                        r and r.HasSubstructMatch(thioether_pattern) for r in reactant_mols if r
                    )

                    if not reactants_have_thioether:
                        thioether_formed = True
                        thioether_depth = depth
                        print(f"Thioether formation detected at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Consider it late-stage if it occurs in the first half of the synthesis
    is_late_stage = (
        thioether_formed and thioether_depth is not None and thioether_depth <= max_depth / 2
    )

    if is_late_stage:
        print(
            f"Late-stage thioether formation detected at depth {thioether_depth} (max depth: {max_depth})"
        )

    return is_late_stage
