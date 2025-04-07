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
    Detects if the synthesis route involves a sequence of sulfur oxidation and reduction
    (thioether → sulfonyl → thioether or similar).
    """
    # Track sulfur transformations
    transformations = []

    def dfs_traverse(node):
        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            try:
                reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles]
                product = Chem.MolFromSmiles(product_smiles)

                if product and all(r for r in reactants):
                    # Check for thioether pattern
                    thioether_pattern = Chem.MolFromSmarts("[#6]-[#16]-[#6]")

                    # Check for sulfonyl pattern
                    sulfonyl_pattern = Chem.MolFromSmarts("[#6]-[#16](=[#8])(=[#8])-[#6]")

                    has_thioether_reactant = any(
                        r.HasSubstructMatch(thioether_pattern) for r in reactants if r
                    )
                    has_sulfonyl_reactant = any(
                        r.HasSubstructMatch(sulfonyl_pattern) for r in reactants if r
                    )

                    has_thioether_product = product.HasSubstructMatch(thioether_pattern)
                    has_sulfonyl_product = product.HasSubstructMatch(sulfonyl_pattern)

                    # Record transformation type
                    if has_thioether_reactant and has_sulfonyl_product:
                        transformations.append("thioether_to_sulfonyl")
                        print("Detected: thioether → sulfonyl")
                    elif has_sulfonyl_reactant and has_thioether_product:
                        transformations.append("sulfonyl_to_thioether")
                        print("Detected: sulfonyl → thioether")
            except:
                print("Error processing reaction SMILES for sulfur transformation detection")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # Check if we have both oxidation and reduction in the sequence
    has_oxidation = "thioether_to_sulfonyl" in transformations
    has_reduction = "sulfonyl_to_thioether" in transformations

    return has_oxidation and has_reduction
