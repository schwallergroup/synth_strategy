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
    This function detects if the synthetic route contains a hydroxamic acid formation
    (coupling of carboxylic acid with hydroxylamine).
    """
    hydroxamic_acid_detected = False

    def dfs_traverse(node):
        nonlocal hydroxamic_acid_detected

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for carboxylic acid and hydroxylamine in reactants
            carboxylic_acid_pattern = Chem.MolFromSmarts("[#6](=[O])[O]")
            hydroxylamine_pattern = Chem.MolFromSmarts("[#7]-[#8]")

            # Check for hydroxamic acid in product
            hydroxamic_acid_pattern = Chem.MolFromSmarts("[#6](=[O])[#7]-[#8]")

            has_carboxylic_acid = False
            has_hydroxylamine = False

            for reactant in reactants:
                mol = Chem.MolFromSmiles(reactant)
                if mol:
                    if mol.HasSubstructMatch(carboxylic_acid_pattern):
                        has_carboxylic_acid = True
                    if mol.HasSubstructMatch(hydroxylamine_pattern):
                        has_hydroxylamine = True

            product_mol = Chem.MolFromSmiles(product)
            has_hydroxamic_acid = product_mol and product_mol.HasSubstructMatch(
                hydroxamic_acid_pattern
            )

            if has_carboxylic_acid and has_hydroxylamine and has_hydroxamic_acid:
                print("Hydroxamic acid formation detected")
                hydroxamic_acid_detected = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return hydroxamic_acid_detected
