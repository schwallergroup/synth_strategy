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
    This function detects if the synthesis route includes multiple Boc protection/deprotection steps.
    """
    boc_operations_count = 0

    def dfs_traverse(node):
        nonlocal boc_operations_count

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for Boc group in reactants or products
                boc_pattern = Chem.MolFromSmarts("[#6]-[#6](-[#6])(-[#6])-[#8]-[#6](=[#8])-[#7]")

                # Check if Boc is being added or removed
                reactants_have_boc = any(
                    Chem.MolFromSmiles(r) and Chem.MolFromSmiles(r).HasSubstructMatch(boc_pattern)
                    for r in reactants
                    if Chem.MolFromSmiles(r)
                )
                product_has_boc = Chem.MolFromSmiles(product) and Chem.MolFromSmiles(
                    product
                ).HasSubstructMatch(boc_pattern)

                if (reactants_have_boc and not product_has_boc) or (
                    not reactants_have_boc and product_has_boc
                ):
                    boc_operations_count += 1
                    print(
                        f"Detected Boc protection/deprotection operation. Total count: {boc_operations_count}"
                    )

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    return boc_operations_count >= 2  # Return True if at least 2 Boc operations are detected
