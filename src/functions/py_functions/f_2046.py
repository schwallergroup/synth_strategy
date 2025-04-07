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
    Detects if the synthesis route involves a late-stage amide coupling (depth 0-1)
    between an amine and an acid chloride or similar activated carboxylic acid.
    """
    found_amide_coupling = False

    def dfs_traverse(node, depth=0):
        nonlocal found_amide_coupling

        if node["type"] == "reaction" and depth <= 1:  # Late stage (depth 0-1)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for amine in reactants
                amine_pattern = Chem.MolFromSmarts("[NH2]")
                # Check for acid chloride in reactants
                acid_chloride_pattern = Chem.MolFromSmarts("[C](=O)[Cl]")
                # Check for amide in product
                amide_pattern = Chem.MolFromSmarts("[NH][C](=O)")

                has_amine = any(
                    Chem.MolFromSmiles(r) is not None
                    and Chem.MolFromSmiles(r).HasSubstructMatch(amine_pattern)
                    for r in reactants
                    if r
                )
                has_acid_chloride = any(
                    Chem.MolFromSmiles(r) is not None
                    and Chem.MolFromSmiles(r).HasSubstructMatch(acid_chloride_pattern)
                    for r in reactants
                    if r
                )
                has_amide = Chem.MolFromSmiles(
                    product
                ) is not None and Chem.MolFromSmiles(product).HasSubstructMatch(
                    amide_pattern
                )

                if has_amine and has_acid_chloride and has_amide:
                    found_amide_coupling = True
                    print(f"Found late-stage amide coupling at depth {depth}")

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    return found_amide_coupling
