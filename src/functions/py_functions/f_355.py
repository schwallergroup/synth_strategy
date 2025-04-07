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
    Detects if the synthesis route uses Heck coupling for carbon backbone construction
    """
    heck_coupling_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal heck_coupling_detected

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for Heck coupling pattern: aryl halide + alkene → aryl-alkene
            aryl_halide_pattern = Chem.MolFromSmarts("[c][Br,Cl,I]")
            alkene_pattern = Chem.MolFromSmarts("[C]=[C]")
            aryl_alkene_pattern = Chem.MolFromSmarts("[c][C]=[C]")

            # Check reactants for aryl halide and alkene
            has_aryl_halide = False
            has_alkene = False

            for r in reactants:
                r_mol = Chem.MolFromSmiles(r)
                if r_mol:
                    if r_mol.HasSubstructMatch(aryl_halide_pattern):
                        has_aryl_halide = True
                    if r_mol.HasSubstructMatch(alkene_pattern):
                        has_alkene = True

            # Check product for aryl-alkene
            product_mol = Chem.MolFromSmiles(product)
            has_aryl_alkene = product_mol and product_mol.HasSubstructMatch(
                aryl_alkene_pattern
            )

            if has_aryl_halide and has_alkene and has_aryl_alkene:
                print(f"Heck coupling detected at depth {depth}")
                heck_coupling_detected = True

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return heck_coupling_detected
