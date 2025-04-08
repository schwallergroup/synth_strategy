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
    This function detects if the synthesis involves a sequence of functional group interconversions:
    carboxylic acid → ketone → olefin with nitrile.
    """
    # Define patterns for functional groups
    carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    ketone_pattern = Chem.MolFromSmarts("[#6][CX3](=O)[#6]")
    cyanoethylidene_pattern = Chem.MolFromSmarts("[#6]=C-C#N")

    # Track if we've seen each transformation
    carboxylic_acid_to_ketone = False
    ketone_to_cyanoethylidene = False

    def dfs_traverse(node):
        nonlocal carboxylic_acid_to_ketone, ketone_to_cyanoethylidene

        if node["type"] == "reaction":
            rsmi = node["metadata"].get("rsmi", "")
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            product_mol = Chem.MolFromSmiles(product)

            # Check for carboxylic acid to ketone transformation
            if product_mol and product_mol.HasSubstructMatch(ketone_pattern):
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(carboxylic_acid_pattern):
                        carboxylic_acid_to_ketone = True
                        print(f"Carboxylic acid to ketone transformation detected: {rsmi}")
                        break

            # Check for ketone to cyanoethylidene transformation
            if product_mol and product_mol.HasSubstructMatch(cyanoethylidene_pattern):
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(ketone_pattern):
                        ketone_to_cyanoethylidene = True
                        print(f"Ketone to cyanoethylidene transformation detected: {rsmi}")
                        break

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    # The strategy is detected if both transformations are present
    strategy_detected = carboxylic_acid_to_ketone and ketone_to_cyanoethylidene
    print(f"Functional group interconversion sequence detected: {strategy_detected}")
    return strategy_detected
