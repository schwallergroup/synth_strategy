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
    This function detects Suzuki coupling reactions in the synthesis route.
    Looks for boronic acid coupling with aryl halide.
    """
    suzuki_found = False

    def dfs_traverse(node):
        nonlocal suzuki_found

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for boronic acid pattern in reactants
            boronic_acid_pattern = Chem.MolFromSmarts("[#5;X3]([#8;X2])[#8;X2]")

            # Check for aryl halide pattern in reactants
            aryl_halide_pattern = Chem.MolFromSmarts("c[F,Cl,Br,I]")

            # Check for biaryl formation in product
            biaryl_pattern = Chem.MolFromSmarts("c:c")

            boronic_acid_present = False
            aryl_halide_present = False

            for reactant in reactants:
                mol = Chem.MolFromSmiles(reactant)
                if mol:
                    if mol.HasSubstructMatch(boronic_acid_pattern):
                        boronic_acid_present = True
                        print(f"Found boronic acid in reactant: {reactant}")
                    if mol.HasSubstructMatch(aryl_halide_pattern):
                        aryl_halide_present = True
                        print(f"Found aryl halide in reactant: {reactant}")

            product_mol = Chem.MolFromSmiles(product)
            biaryl_present = False
            if product_mol and product_mol.HasSubstructMatch(biaryl_pattern):
                biaryl_present = True
                print(f"Found biaryl in product: {product}")

            if boronic_acid_present and aryl_halide_present and biaryl_present:
                suzuki_found = True
                print("Suzuki coupling detected!")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return suzuki_found
