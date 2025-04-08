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
    This function detects transformation of aryl halide to boronic acid.
    """
    halogen_to_boronic_detected = False

    def dfs_traverse(node):
        nonlocal halogen_to_boronic_detected

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for aryl halide in reactants
            aryl_halide_present = False
            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if reactant_mol:
                    aryl_halide_pattern = Chem.MolFromSmarts("[c][Br,Cl,I]")
                    if reactant_mol.HasSubstructMatch(aryl_halide_pattern):
                        aryl_halide_present = True
                        break

            # Check for boronic acid in product
            product_mol = Chem.MolFromSmiles(product)
            if product_mol and aryl_halide_present:
                boronic_acid_pattern = Chem.MolFromSmarts("[c][B]([OH])[OH]")
                if product_mol.HasSubstructMatch(boronic_acid_pattern):
                    halogen_to_boronic_detected = True
                    print("Halogen to boronic acid transformation detected in reaction:", rsmi)

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return halogen_to_boronic_detected
