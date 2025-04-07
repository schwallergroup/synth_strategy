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
    This function detects if the synthesis includes a late-stage nitro reduction to amine.
    Late stage is defined as occurring in the first half of the synthesis (depth < max_depth/2).
    """
    found_nitro_reduction = False
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal found_nitro_reduction, max_depth

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this is a nitro reduction reaction
            nitro_pattern = Chem.MolFromSmarts("[N+](=O)[O-]")
            amine_pattern = Chem.MolFromSmarts("[NH2]")

            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]
            product_mol = Chem.MolFromSmiles(product)

            if product_mol and all(r for r in reactant_mols):
                has_nitro_in_reactants = any(
                    r.HasSubstructMatch(nitro_pattern) for r in reactant_mols if r
                )
                has_amine_in_product = (
                    product_mol.HasSubstructMatch(amine_pattern) if product_mol else False
                )

                if (
                    has_nitro_in_reactants and has_amine_in_product and depth < 2
                ):  # Late stage (depth < 2)
                    found_nitro_reduction = True
                    print(f"Found late-stage nitro reduction at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return found_nitro_reduction
