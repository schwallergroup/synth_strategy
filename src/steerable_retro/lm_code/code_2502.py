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
    This function detects a synthetic strategy involving conversion of aryl halides
    to boronic acid/ester derivatives.
    """
    found_halogen_to_boryl = False

    def dfs_traverse(node):
        nonlocal found_halogen_to_boryl

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for aryl halide in reactants
            aryl_halide_pattern = Chem.MolFromSmarts("[#6]:[#6]-[#9,#17,#35,#53]")

            # Check for boronic acid/ester in product
            boronic_pattern = Chem.MolFromSmarts("[#6]-[#5](-[#8])-[#8]")

            try:
                # Check if any reactant has aryl halide
                has_aryl_halide = False
                for reactant in reactants:
                    react_mol = Chem.MolFromSmiles(reactant)
                    if react_mol and react_mol.HasSubstructMatch(aryl_halide_pattern):
                        has_aryl_halide = True
                        break

                # Check if product has boronic acid/ester
                prod_mol = Chem.MolFromSmiles(product)
                has_boronic = prod_mol and prod_mol.HasSubstructMatch(boronic_pattern)

                if has_aryl_halide and has_boronic:
                    found_halogen_to_boryl = True
                    print(f"Detected halogen to boryl transformation: {rsmi}")
            except:
                pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return found_halogen_to_boryl
