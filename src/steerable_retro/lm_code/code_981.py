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
    Detects if the synthesis route includes a transformation from
    carboxylic acid to amide.
    """
    transformation_detected = False

    def dfs_traverse(node):
        nonlocal transformation_detected

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for carboxylic acid in reactants
                carboxylic_acid_pattern = Chem.MolFromSmarts("[C](=O)[OH]")

                # Check for amide in product
                amide_pattern = Chem.MolFromSmarts("[C](=O)[NH2]")

                # Check if transformation occurred
                if any(
                    Chem.MolFromSmiles(r)
                    and Chem.MolFromSmiles(r).HasSubstructMatch(carboxylic_acid_pattern)
                    for r in reactants
                    if Chem.MolFromSmiles(r)
                ):
                    prod_mol = Chem.MolFromSmiles(product)
                    if prod_mol and prod_mol.HasSubstructMatch(amide_pattern):
                        print("Carboxylic acid to amide transformation detected")
                        transformation_detected = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return transformation_detected
