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
    Detects if the synthesis route involves multiple nitrogen deprotection steps
    (e.g., Boc deprotection, sulfonyl deprotection).
    """
    deprotection_count = 0

    def dfs_traverse(node):
        nonlocal deprotection_count

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_str = rsmi.split(">")[0]
                product_str = rsmi.split(">")[-1]

                # Boc deprotection pattern
                boc_protected = Chem.MolFromSmarts("[NH][C](=[O])[O][C]([C])([C])[C]")
                free_amine = Chem.MolFromSmarts("[NH2]")

                # Sulfonyl deprotection pattern
                sulfonyl_protected = Chem.MolFromSmarts("[n][S](=[O])(=[O])[c]")
                free_indole_nh = Chem.MolFromSmarts("[nH]")

                reactants = [Chem.MolFromSmiles(r) for r in reactants_str.split(".")]
                product = Chem.MolFromSmiles(product_str)

                # Check for Boc deprotection
                has_boc_reactant = any(
                    r is not None and r.HasSubstructMatch(boc_protected)
                    for r in reactants
                )
                has_free_amine_product = (
                    product is not None and product.HasSubstructMatch(free_amine)
                )

                # Check for sulfonyl deprotection
                has_sulfonyl_reactant = any(
                    r is not None and r.HasSubstructMatch(sulfonyl_protected)
                    for r in reactants
                )
                has_free_indole_product = (
                    product is not None and product.HasSubstructMatch(free_indole_nh)
                )

                if (has_boc_reactant and has_free_amine_product) or (
                    has_sulfonyl_reactant and has_free_indole_product
                ):
                    print("Found nitrogen deprotection step")
                    deprotection_count += 1

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return deprotection_count >= 2  # At least two deprotection steps
