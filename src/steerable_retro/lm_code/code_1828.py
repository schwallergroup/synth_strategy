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
    Detects if the synthetic route preserves certain structural motifs (chlorobenzene,
    isoxazole, and morpholine) throughout the synthesis.
    """
    # Define the motifs to track
    motifs = {
        "chlorobenzene": Chem.MolFromSmarts("[#6]1:[#6]:[#6](-[#17]):[#6]:[#6]:[#6]:1"),
        "isoxazole": Chem.MolFromSmarts("[#6]1:[#7]:[#8]:[#6]:[#6]:1"),
        "morpholine": Chem.MolFromSmarts("[#7]1-[#6]-[#6]-[#8]-[#6]-[#6]-1"),
    }

    # Track which motifs are present at each step
    motifs_present = {name: [] for name in motifs}
    all_products = []

    def dfs_traverse(node):
        if node["type"] == "reaction":
            # Extract product
            rsmi = node["metadata"]["rsmi"]
            product_smiles = rsmi.split(">")[-1]
            product = Chem.MolFromSmiles(product_smiles) if product_smiles else None

            if product:
                all_products.append(product)
                # Check each motif
                for name, pattern in motifs.items():
                    if product.HasSubstructMatch(pattern):
                        motifs_present[name].append(True)
                    else:
                        motifs_present[name].append(False)

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Check if we have products and if motifs are preserved
    if not all_products:
        return False

    preserved_motifs = []
    for name, presence_list in motifs_present.items():
        if presence_list and all(presence_list):
            preserved_motifs.append(name)
            print(f"Motif {name} is preserved throughout synthesis")

    # Strategy is present if at least two motifs are preserved
    strategy_present = len(preserved_motifs) >= 2
    print(f"Preserved motifs strategy detected: {strategy_present}")
    return strategy_present
