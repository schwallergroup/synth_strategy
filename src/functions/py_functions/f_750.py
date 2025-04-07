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
    Detects a synthetic strategy involving Williamson ether synthesis (phenol + alkyl halide).
    """
    has_williamson_ether = False

    def dfs_traverse(node):
        nonlocal has_williamson_ether

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                if reactants and product:
                    # Check for phenol and alkyl halide patterns in reactants
                    phenol_pattern = Chem.MolFromSmarts("[OH][c]")
                    alkyl_halide_pattern = Chem.MolFromSmarts("[Br,Cl,I][CH2][CH2]")

                    has_phenol = False
                    has_alkyl_halide = False

                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol:
                            if reactant_mol.HasSubstructMatch(phenol_pattern):
                                has_phenol = True
                            if reactant_mol.HasSubstructMatch(alkyl_halide_pattern):
                                has_alkyl_halide = True

                    # Check for aryl ether pattern in product
                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol:
                        aryl_ether_pattern = Chem.MolFromSmarts("[c][O][CH2][CH2]")
                        if product_mol.HasSubstructMatch(aryl_ether_pattern):
                            if has_phenol and has_alkyl_halide:
                                has_williamson_ether = True

        # Recursively process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Williamson ether synthesis strategy: {has_williamson_ether}")
    return has_williamson_ether
