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
    Detects if the synthesis route involves multiple interconversions between
    carboxylic acid derivatives (ester, acid, amide).
    """
    interconversion_count = 0

    def dfs_traverse(node):
        nonlocal interconversion_count

        if node["type"] == "reaction" and node.get("metadata", {}).get("rsmi"):
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0]
            product_smiles = rsmi.split(">")[-1]

            # Define patterns for carboxylic acid derivatives
            acid_pattern = Chem.MolFromSmarts("[C](=[O])[OH]")
            ester_pattern = Chem.MolFromSmarts("[C](=[O])[O][C]")
            amide_pattern = Chem.MolFromSmarts("[C](=[O])[N]")

            reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles.split(".")]
            product = Chem.MolFromSmiles(product_smiles)

            # Check for acid in reactants
            acid_in_reactants = any(
                r is not None and r.HasSubstructMatch(acid_pattern) for r in reactants
            )

            # Check for ester in reactants
            ester_in_reactants = any(
                r is not None and r.HasSubstructMatch(ester_pattern) for r in reactants
            )

            # Check for amide in reactants
            amide_in_reactants = any(
                r is not None and r.HasSubstructMatch(amide_pattern) for r in reactants
            )

            # Check for acid in product
            acid_in_product = product is not None and product.HasSubstructMatch(
                acid_pattern
            )

            # Check for ester in product
            ester_in_product = product is not None and product.HasSubstructMatch(
                ester_pattern
            )

            # Check for amide in product
            amide_in_product = product is not None and product.HasSubstructMatch(
                amide_pattern
            )

            # Check for interconversions
            if (
                (acid_in_reactants and (ester_in_product or amide_in_product))
                or (ester_in_reactants and (acid_in_product or amide_in_product))
                or (amide_in_reactants and (acid_in_product or ester_in_product))
            ):
                interconversion_count += 1
                print(
                    f"Detected carboxylic acid derivative interconversion in reaction: {node.get('metadata', {}).get('ID', '')}"
                )

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return interconversion_count >= 2  # At least 2 interconversions
