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
    This function detects synthetic routes involving the conversion of terminal alcohols
    to halides (particularly chlorides, bromides, or iodides).
    """
    # Track if we found evidence of alcohol to halide conversion
    has_conversion = False

    def dfs_traverse(node):
        nonlocal has_conversion

        if node["type"] == "reaction":
            # Extract reactants and products
            rsmi = node["metadata"].get("rsmi", "")
            if not rsmi:
                return

            parts = rsmi.split(">")
            if len(parts) < 3:
                return

            reactants = parts[0].split(".")
            products = parts[2].split(".")

            # Check for alcohol to halide conversion
            for reactant in reactants:
                try:
                    r_mol = Chem.MolFromSmiles(reactant)
                    if r_mol:
                        # Check for terminal alcohol
                        terminal_alcohol = Chem.MolFromSmarts("[#6]-[#8;H1]")
                        if r_mol.HasSubstructMatch(terminal_alcohol):
                            # Look for corresponding product with terminal halide
                            for product in products:
                                try:
                                    p_mol = Chem.MolFromSmiles(product)
                                    if p_mol:
                                        # Check for terminal halide (Cl, Br, I)
                                        terminal_halide = Chem.MolFromSmarts("[#6]-[#17,#35,#53]")
                                        if p_mol.HasSubstructMatch(terminal_halide):
                                            has_conversion = True
                                            print(
                                                f"Found alcohol to halide conversion: {reactant} -> {product}"
                                            )
                                except:
                                    continue
                except:
                    continue

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Alcohol to halide conversion strategy detected: {has_conversion}")
    return has_conversion
