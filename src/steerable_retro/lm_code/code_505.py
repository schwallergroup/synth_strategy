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
    This function detects if the synthetic route involves multiple interconversions
    of carboxylic acid derivatives (acid, ester, amide).
    """
    transformations = []

    def dfs_traverse(node):
        nonlocal transformations

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for carboxylic acid derivative transformations
            acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
            ester_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2][C]")
            amide_pattern = Chem.MolFromSmarts("[CX3](=O)[NX3]")

            reactant_has_acid = any(
                Chem.MolFromSmiles(r).HasSubstructMatch(acid_pattern) for r in reactants
            )
            reactant_has_ester = any(
                Chem.MolFromSmiles(r).HasSubstructMatch(ester_pattern) for r in reactants
            )
            reactant_has_amide = any(
                Chem.MolFromSmiles(r).HasSubstructMatch(amide_pattern) for r in reactants
            )

            product_has_acid = Chem.MolFromSmiles(product).HasSubstructMatch(acid_pattern)
            product_has_ester = Chem.MolFromSmiles(product).HasSubstructMatch(ester_pattern)
            product_has_amide = Chem.MolFromSmiles(product).HasSubstructMatch(amide_pattern)

            # Detect transformations
            if reactant_has_acid and product_has_ester:
                transformations.append("acid_to_ester")
                print("Transformation: acid to ester")
            elif reactant_has_acid and product_has_amide:
                transformations.append("acid_to_amide")
                print("Transformation: acid to amide")
            elif reactant_has_ester and product_has_acid:
                transformations.append("ester_to_acid")
                print("Transformation: ester to acid")
            elif reactant_has_ester and product_has_amide:
                transformations.append("ester_to_amide")
                print("Transformation: ester to amide")
            elif reactant_has_amide and product_has_acid:
                transformations.append("amide_to_acid")
                print("Transformation: amide to acid")
            elif reactant_has_amide and product_has_ester:
                transformations.append("amide_to_ester")
                print("Transformation: amide to ester")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if multiple carboxylic acid derivative transformations are found
    return len(transformations) >= 2
