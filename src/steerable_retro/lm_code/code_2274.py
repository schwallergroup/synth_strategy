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
    This function detects a sequence of functional group interconversions including
    ester, carboxylic acid, alcohol, and aldehyde.
    """
    # Track functional group transformations
    transformations = []

    def dfs_traverse(node):
        nonlocal transformations

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0]
                product = rsmi.split(">")[-1]

                reactant_mol = Chem.MolFromSmiles(reactants)
                product_mol = Chem.MolFromSmiles(product)

                if reactant_mol and product_mol:
                    # Define patterns for functional groups
                    ester_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2][CX4]")
                    acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
                    alcohol_pattern = Chem.MolFromSmarts("[CX4][OX2H]")
                    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)")

                    # Check for transformations
                    if reactant_mol.HasSubstructMatch(
                        ester_pattern
                    ) and product_mol.HasSubstructMatch(acid_pattern):
                        transformations.append("ester_to_acid")
                        print("Detected ester to acid transformation")

                    if reactant_mol.HasSubstructMatch(
                        acid_pattern
                    ) and product_mol.HasSubstructMatch(alcohol_pattern):
                        transformations.append("acid_to_alcohol")
                        print("Detected acid to alcohol transformation")

                    if reactant_mol.HasSubstructMatch(
                        alcohol_pattern
                    ) and product_mol.HasSubstructMatch(aldehyde_pattern):
                        transformations.append("alcohol_to_aldehyde")
                        print("Detected alcohol to aldehyde transformation")

                    if reactant_mol.HasSubstructMatch(
                        aldehyde_pattern
                    ) and product_mol.HasSubstructMatch(acid_pattern):
                        transformations.append("aldehyde_to_acid")
                        print("Detected aldehyde to acid transformation")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Check if we have at least 3 different transformations
    unique_transformations = set(transformations)
    return len(unique_transformations) >= 3
