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
    Detects if the synthesis route has a late-stage amide coupling reaction
    (depth 0 or 1) that forms a key C-N bond.
    """
    found_amide_coupling = False

    def dfs_traverse(node, depth=0):
        nonlocal found_amide_coupling

        if node["type"] == "reaction" and depth <= 1:  # Late stage (depth 0 or 1)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_str = rsmi.split(">")[0]
                product_str = rsmi.split(">")[-1]

                # Check for acid chloride or similar activated carboxylic acid
                acid_pattern = Chem.MolFromSmarts("[C](=[O])[Cl,Br,I,O]")
                # Check for amine nucleophile
                amine_pattern = Chem.MolFromSmarts("[N;H1,H2]")
                # Check for amide in product
                amide_pattern = Chem.MolFromSmarts("[C](=[O])[N]")

                reactants = [Chem.MolFromSmiles(r) for r in reactants_str.split(".")]
                product = Chem.MolFromSmiles(product_str)

                # Check if reactants contain acid chloride and amine
                has_acid = any(
                    r is not None and r.HasSubstructMatch(acid_pattern) for r in reactants
                )
                has_amine = any(
                    r is not None and r.HasSubstructMatch(amine_pattern) for r in reactants
                )

                # Check if product contains amide
                has_amide_product = product is not None and product.HasSubstructMatch(amide_pattern)

                if has_acid and has_amine and has_amide_product:
                    print(f"Found late-stage amide coupling at depth {depth}")
                    found_amide_coupling = True

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return found_amide_coupling
