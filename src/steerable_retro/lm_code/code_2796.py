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
    Detects if the synthesis follows a linear strategy with carbamate formation.
    """
    # Track if we found the pattern
    found_pattern = False
    # Track the depth at which we find key transformations
    carbamate_formation_depth = None
    fragment_counts = []

    def dfs_traverse(node, depth=0):
        nonlocal found_pattern, carbamate_formation_depth, fragment_counts

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Count fragments
            fragment_counts.append(len(reactants))

            # Check for carbamate formation
            carbamate_pattern = Chem.MolFromSmarts("[#8][#6](=[#8])[#7]")

            product_mol = Chem.MolFromSmiles(product)
            if product_mol and product_mol.HasSubstructMatch(carbamate_pattern):
                # Check if any reactant has chloroformate
                chloroformate_pattern = Chem.MolFromSmarts("[Cl][#6](=[#8])[#8]")
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(chloroformate_pattern):
                        carbamate_formation_depth = depth
                        print(f"Found carbamate formation at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if we found the pattern
    if carbamate_formation_depth is not None:
        # Check if synthesis is linear (no more than 2 fragments per step)
        if all(count <= 2 for count in fragment_counts):
            found_pattern = True
            print("Found linear synthesis with carbamate formation")

    return found_pattern
