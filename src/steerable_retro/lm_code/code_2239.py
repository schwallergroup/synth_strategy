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
    This function detects if the synthesis uses an early-stage Sonogashira coupling
    for C-C bond formation via alkyne coupling.
    """
    sonogashira_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal sonogashira_detected

        if node["type"] == "reaction" and depth >= 3:  # Early in synthesis (high depth)
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for patterns indicative of Sonogashira coupling
            halogen_aromatic_pattern = Chem.MolFromSmarts("[Br,I][c,n]")
            terminal_alkyne_pattern = Chem.MolFromSmarts("[C]#[CH]")
            coupled_alkyne_pattern = Chem.MolFromSmarts("[c,n][C]#[C][c,n]")

            # Check reactants for halogen-aromatic and terminal alkyne
            has_halogen_aromatic = False
            has_terminal_alkyne = False

            for reactant in reactants:
                mol = Chem.MolFromSmiles(reactant)
                if mol:
                    if mol.HasSubstructMatch(halogen_aromatic_pattern):
                        has_halogen_aromatic = True
                    if mol.HasSubstructMatch(terminal_alkyne_pattern):
                        has_terminal_alkyne = True

            # Check product for coupled alkyne
            product_mol = Chem.MolFromSmiles(product)
            has_coupled_alkyne = product_mol and product_mol.HasSubstructMatch(
                coupled_alkyne_pattern
            )

            if has_halogen_aromatic and has_terminal_alkyne and has_coupled_alkyne:
                sonogashira_detected = True
                print(f"Early-stage Sonogashira coupling detected at depth {depth}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    return sonogashira_detected
