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
    This function detects a linear synthesis strategy focused on sequential
    modifications at a benzylic position while maintaining the aromatic core.
    """
    # Track benzylic modifications
    benzylic_modifications = 0
    total_reactions = 0

    def dfs_traverse(node):
        nonlocal benzylic_modifications, total_reactions

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            total_reactions += 1

            # Check for benzylic position modifications
            benzylic_patterns = [
                Chem.MolFromSmarts("[c][C](=O)[OH]"),  # Benzylic carboxylic acid
                Chem.MolFromSmarts("[c][C][OH]"),  # Benzylic alcohol
                Chem.MolFromSmarts("[c][C][Cl]"),  # Benzylic chloride
                Chem.MolFromSmarts("[c][C][O][C]"),  # Benzylic ether
            ]

            product_mol = Chem.MolFromSmiles(product)
            if product_mol:
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if Chem.MolFromSmiles(r)]

                # Check if any reactant has a different benzylic group than the product
                reactant_matches = set()
                for mol in reactant_mols:
                    if mol:
                        for pattern in benzylic_patterns:
                            if mol.HasSubstructMatch(pattern):
                                reactant_matches.add(pattern)

                product_matches = set()
                for pattern in benzylic_patterns:
                    if product_mol.HasSubstructMatch(pattern):
                        product_matches.add(pattern)

                # If the benzylic group changed, count it as a modification
                if reactant_matches != product_matches and (reactant_matches or product_matches):
                    benzylic_modifications += 1
                    print(f"Found benzylic modification in reaction {total_reactions}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Strategy is present if we have at least 2 benzylic modifications in a linear synthesis
    result = benzylic_modifications >= 2 and total_reactions >= 3
    print(
        f"Linear benzylic modification strategy detected: {result} ({benzylic_modifications}/{total_reactions} reactions)"
    )
    return result
