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
    Detects Suzuki coupling for biaryl formation in the synthetic route.
    Looks for boronic acid and aryl halide reactants with biaryl product.
    """
    suzuki_detected = False

    def dfs_traverse(node):
        nonlocal suzuki_detected

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for boronic acid pattern
            boronic_acid_pattern = Chem.MolFromSmarts("[OX2H]-[BX3](-[OX2H])-[c,C]")

            # Check for aryl halide pattern (I, Br preferred for Suzuki)
            aryl_halide_pattern = Chem.MolFromSmarts("[c]-[I,Br,Cl]")

            # Check for biaryl pattern in product
            biaryl_pattern = Chem.MolFromSmarts("c:c-c:c")

            has_boronic_acid = False
            has_aryl_halide = False

            for reactant in reactants:
                mol = Chem.MolFromSmiles(reactant)
                if mol:
                    if mol.HasSubstructMatch(boronic_acid_pattern):
                        has_boronic_acid = True
                    if mol.HasSubstructMatch(aryl_halide_pattern):
                        has_aryl_halide = True

            product_mol = Chem.MolFromSmiles(product)
            has_biaryl = False
            if product_mol and product_mol.HasSubstructMatch(biaryl_pattern):
                has_biaryl = True

            if has_boronic_acid and has_aryl_halide and has_biaryl:
                print("Suzuki coupling detected for biaryl formation")
                suzuki_detected = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return suzuki_detected
