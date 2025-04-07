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
    This function detects benzyl ether formation as a linking strategy.
    """
    found = False

    def dfs_traverse(node, depth=0):
        nonlocal found

        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if product contains benzyl ether
                product_mol = Chem.MolFromSmiles(product)
                if product_mol and product_mol.HasSubstructMatch(
                    Chem.MolFromSmarts("[c][CH2][O][c]")
                ):

                    # Check if one reactant is benzyl halide
                    benzyl_halide_pattern = Chem.MolFromSmarts("[c][CH2][Cl,Br,I]")

                    # Check if one reactant is phenol
                    phenol_pattern = Chem.MolFromSmarts("[c][OH]")

                    has_benzyl_halide = False
                    has_phenol = False

                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol:
                            if reactant_mol.HasSubstructMatch(benzyl_halide_pattern):
                                has_benzyl_halide = True
                            if reactant_mol.HasSubstructMatch(phenol_pattern):
                                has_phenol = True

                    if has_benzyl_halide and has_phenol:
                        found = True
                        print(f"Found benzyl ether formation at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return found
