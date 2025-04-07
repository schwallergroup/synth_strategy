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
    Detects if the synthesis route employs a late-stage Suzuki coupling (depth 0-1)
    to connect two aromatic fragments.
    """
    suzuki_detected = False
    max_depth = 1  # Consider only reactions at depth 0-1 as "late-stage"

    def dfs_traverse(node, depth=0):
        nonlocal suzuki_detected

        if node["type"] == "reaction" and depth <= max_depth:
            # Extract reactants and product
            rsmi = node["metadata"].get("rsmi", "")
            if not rsmi:
                return

            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for Suzuki coupling patterns
            boronic_acid_pattern = Chem.MolFromSmarts("[c]-[B]([O])[O]")
            aryl_halide_pattern = Chem.MolFromSmarts("[c]-[#17,#35,#53]")
            biaryl_pattern = Chem.MolFromSmarts("[c]-[c]")

            # Check if reactants contain boronic acid and aryl halide
            has_boronic_acid = False
            has_aryl_halide = False

            for reactant in reactants:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(boronic_acid_pattern):
                        has_boronic_acid = True
                    if mol and mol.HasSubstructMatch(aryl_halide_pattern):
                        has_aryl_halide = True
                except:
                    continue

            # Check if product has biaryl bond
            try:
                prod_mol = Chem.MolFromSmiles(product)
                has_biaryl = prod_mol and prod_mol.HasSubstructMatch(biaryl_pattern)
            except:
                has_biaryl = False

            # If all conditions are met, it's likely a Suzuki coupling
            if has_boronic_acid and has_aryl_halide and has_biaryl:
                print(f"Detected Suzuki coupling at depth {depth}")
                suzuki_detected = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return suzuki_detected
