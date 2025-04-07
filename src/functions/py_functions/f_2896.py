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
    Detects if the synthesis involves oxidation of a thioether to a sulfoxide.
    """
    sulfur_oxidation_found = False

    def dfs_traverse(node, depth=0):
        nonlocal sulfur_oxidation_found

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for sulfur oxidation
            product_mol = Chem.MolFromSmiles(product)
            reactant_mols = [
                Chem.MolFromSmiles(r) for r in reactants if Chem.MolFromSmiles(r)
            ]

            if not product_mol or not reactant_mols:
                return

            # Check if product has sulfoxide but reactant has thioether
            sulfoxide_pattern = Chem.MolFromSmarts("[#16](=[#8])")
            thioether_pattern = Chem.MolFromSmarts("[#6]-[#16]-[#6]")

            if product_mol.HasSubstructMatch(sulfoxide_pattern):
                for reactant_mol in reactant_mols:
                    if reactant_mol.HasSubstructMatch(
                        thioether_pattern
                    ) and not reactant_mol.HasSubstructMatch(sulfoxide_pattern):
                        sulfur_oxidation_found = True
                        print(f"Found sulfur oxidation at depth {depth}")
                        break

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return sulfur_oxidation_found
