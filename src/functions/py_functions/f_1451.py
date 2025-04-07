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
    This function detects if the synthetic route involves a carboxylic acid → acid chloride → amide
    functional group transformation sequence.
    """
    # Flags to track if we found each step in the sequence
    found_acid_to_acid_chloride = False
    found_acid_chloride_to_amide = False

    def dfs_traverse(node):
        nonlocal found_acid_to_acid_chloride, found_acid_chloride_to_amide

        # Only process reaction nodes
        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            # Get reaction SMILES
            rsmi = node["metadata"]["rsmi"]

            # Extract reactants and products
            reactants_str = rsmi.split(">")[0]
            product_str = rsmi.split(">")[-1]

            # Parse reactants and product
            reactants = [Chem.MolFromSmiles(r) for r in reactants_str.split(".") if r]
            product = Chem.MolFromSmiles(product_str)

            if product and all(r for r in reactants):
                # Check for carboxylic acid pattern in reactants
                acid_pattern = Chem.MolFromSmarts("[#6](=O)[OH]")
                has_acid = any(
                    r.HasSubstructMatch(acid_pattern) for r in reactants if r
                )

                # Check for acid chloride pattern in product
                acid_chloride_pattern = Chem.MolFromSmarts("[#6](=O)[Cl]")
                has_acid_chloride_product = (
                    product.HasSubstructMatch(acid_chloride_pattern)
                    if product
                    else False
                )

                # Check for acid chloride pattern in reactants
                has_acid_chloride_reactant = any(
                    r.HasSubstructMatch(acid_chloride_pattern) for r in reactants if r
                )

                # Check for amine pattern in reactants
                amine_pattern = Chem.MolFromSmarts("[NH2]")
                has_amine = any(
                    r.HasSubstructMatch(amine_pattern) for r in reactants if r
                )

                # Check for amide pattern in product
                amide_pattern = Chem.MolFromSmarts("[#6](=O)[#7]")
                has_amide_product = (
                    product.HasSubstructMatch(amide_pattern) if product else False
                )

                # Check for acid to acid chloride conversion
                if has_acid and has_acid_chloride_product:
                    print("Found carboxylic acid to acid chloride conversion")
                    found_acid_to_acid_chloride = True

                # Check for acid chloride to amide conversion
                if has_acid_chloride_reactant and has_amine and has_amide_product:
                    print("Found acid chloride to amide conversion")
                    found_acid_chloride_to_amide = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True only if we found the complete sequence
    return found_acid_to_acid_chloride and found_acid_chloride_to_amide
