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
    This function detects if the synthetic route preserves a methoxy group
    on an aromatic ring throughout the synthesis.
    """
    methoxy_reactions_count = 0
    total_reactions = 0

    def dfs_traverse(node):
        nonlocal methoxy_reactions_count, total_reactions

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                total_reactions += 1
                rsmi = node["metadata"]["rsmi"]
                product = rsmi.split(">")[-1]

                # Check for methoxy-aromatic pattern in product
                methoxy_aromatic_pattern = Chem.MolFromSmarts("c[OD2][CH3]")

                try:
                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol and product_mol.HasSubstructMatch(methoxy_aromatic_pattern):
                        print(f"Found methoxy-aromatic in reaction product: {product}")
                        methoxy_reactions_count += 1
                except:
                    print(f"Error processing SMILES: {product}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # Return True if methoxy-aromatic is present in at least 80% of reactions
    return total_reactions > 0 and (methoxy_reactions_count / total_reactions) >= 0.8
