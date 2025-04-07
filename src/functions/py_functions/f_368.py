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
    This function detects if the synthesis route involves a late-stage amide coupling
    (depth 0 or 1) between an acid chloride and an amine.
    """
    amide_formation_found = False

    def dfs_traverse(node, depth=0):
        nonlocal amide_formation_found

        if node["type"] == "reaction" and depth <= 1:  # Only check late-stage reactions
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for acid chloride in reactants
                acid_chloride_pattern = Chem.MolFromSmarts("[C](=O)Cl")
                # Check for primary amine in reactants
                amine_pattern = Chem.MolFromSmarts("[NH2]")
                # Check for amide in product
                amide_pattern = Chem.MolFromSmarts("[NH][C](=O)")

                acid_chloride_present = any(
                    Chem.MolFromSmiles(r) is not None
                    and Chem.MolFromSmiles(r).HasSubstructMatch(acid_chloride_pattern)
                    for r in reactants
                    if r
                )
                amine_present = any(
                    Chem.MolFromSmiles(r) is not None
                    and Chem.MolFromSmiles(r).HasSubstructMatch(amine_pattern)
                    for r in reactants
                    if r
                )
                amide_in_product = Chem.MolFromSmiles(
                    product
                ) is not None and Chem.MolFromSmiles(product).HasSubstructMatch(
                    amide_pattern
                )

                if acid_chloride_present and amine_present and amide_in_product:
                    print(f"Found amide formation at depth {depth}")
                    amide_formation_found = True

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    return amide_formation_found
