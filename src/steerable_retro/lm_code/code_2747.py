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
    Detects a carboxylic acid protection-deprotection strategy
    """
    has_protection = False
    has_deprotection = False

    def dfs_traverse(node):
        nonlocal has_protection, has_deprotection

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
            product = Chem.MolFromSmiles(product_smiles) if product_smiles else None

            if product and reactants:
                # Check for carboxylic acid protection (acid to ester)
                acid_pattern = Chem.MolFromSmarts("[C](=[O])[OH]")
                ester_pattern = Chem.MolFromSmarts("[C](=[O])[O][C]")

                if any(
                    r and r.HasSubstructMatch(acid_pattern) for r in reactants
                ) and product.HasSubstructMatch(ester_pattern):
                    print("Detected carboxylic acid protection")
                    has_protection = True

                # Check for carboxylic acid deprotection (ester to acid)
                if any(
                    r and r.HasSubstructMatch(ester_pattern) for r in reactants
                ) and product.HasSubstructMatch(acid_pattern):
                    print("Detected carboxylic acid deprotection")
                    has_deprotection = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Strategy is present if we have both protection and deprotection
    return has_protection and has_deprotection
