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
    This function detects a synthetic strategy involving TMS protection and deprotection
    of terminal alkynes in the synthesis route.
    """
    has_tms_protection = False
    has_tms_deprotection = False

    def dfs_traverse(node):
        nonlocal has_tms_protection, has_tms_deprotection

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for TMS protection (terminal alkyne → TMS-alkyne)
            terminal_alkyne_pattern = Chem.MolFromSmarts("[C]#[CH]")
            tms_alkyne_pattern = Chem.MolFromSmarts("[C]#[C][Si](C)(C)C")

            # Check for TMS protection
            has_terminal_alkyne_reactant = False
            for reactant in reactants:
                mol = Chem.MolFromSmiles(reactant)
                if mol and mol.HasSubstructMatch(terminal_alkyne_pattern):
                    has_terminal_alkyne_reactant = True
                    break

            product_mol = Chem.MolFromSmiles(product)
            if (
                has_terminal_alkyne_reactant
                and product_mol
                and product_mol.HasSubstructMatch(tms_alkyne_pattern)
            ):
                has_tms_protection = True
                print(f"Found TMS protection reaction: {rsmi}")

            # Check for TMS deprotection (TMS-alkyne → terminal alkyne)
            has_tms_alkyne_reactant = False
            for reactant in reactants:
                mol = Chem.MolFromSmiles(reactant)
                if mol and mol.HasSubstructMatch(tms_alkyne_pattern):
                    has_tms_alkyne_reactant = True
                    break

            if (
                has_tms_alkyne_reactant
                and product_mol
                and product_mol.HasSubstructMatch(terminal_alkyne_pattern)
            ):
                has_tms_deprotection = True
                print(f"Found TMS deprotection reaction: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if both TMS protection and deprotection are found
    result = has_tms_protection and has_tms_deprotection
    print(f"TMS alkyne protection-deprotection strategy detected: {result}")
    return result
