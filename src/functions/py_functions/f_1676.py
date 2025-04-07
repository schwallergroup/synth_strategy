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
    Detects amide formation reactions in the synthesis route.
    Looks for reactions where an amine and carboxylic acid/derivative form an amide.
    """
    result = False

    def dfs_traverse(node, depth=0):
        nonlocal result

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0]
            product_smiles = rsmi.split(">")[-1]

            # Parse reactants and product
            reactants = [
                Chem.MolFromSmiles(r) for r in reactants_smiles.split(".") if r
            ]
            product = Chem.MolFromSmiles(product_smiles)

            if not all(reactants) or not product:
                print("Warning: Could not parse all molecules in reaction")
                return

            # Check for amine in reactants
            amine_pattern = Chem.MolFromSmarts("[#6]-[NH2]")
            has_amine = any(mol.HasSubstructMatch(amine_pattern) for mol in reactants)

            # Check for carboxylic acid or derivative in reactants
            acid_pattern = Chem.MolFromSmarts("[#6]-C(=O)-[OH]")
            ester_pattern = Chem.MolFromSmarts("[#6]-C(=O)-O-[#6]")
            has_acid_or_ester = any(
                mol.HasSubstructMatch(acid_pattern)
                or mol.HasSubstructMatch(ester_pattern)
                for mol in reactants
            )

            # Check for amide in product
            amide_pattern = Chem.MolFromSmarts("[#6]-C(=O)-N")
            has_amide_product = product.HasSubstructMatch(amide_pattern)

            if has_amine and has_acid_or_ester and has_amide_product:
                print(f"Found amide formation at depth {depth}")
                result = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return result
