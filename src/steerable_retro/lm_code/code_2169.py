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
    Detects a strategy where an alcohol is first protected, then another alcohol
    is converted to a leaving group (mesylate) and then to an iodide.
    """
    # Initialize tracking variables
    alcohol_protection_steps = 0
    alcohol_to_leaving_group = False
    leaving_group_substitution = False

    # SMARTS patterns
    tbdms_pattern = Chem.MolFromSmarts("[OX2][Si]([C])([C])[C]([C])([C])[C]")
    benzyl_pattern = Chem.MolFromSmarts("[OX2][CH2][cX3]1[cX3][cX3][cX3][cX3][cX3]1")
    mesylate_pattern = Chem.MolFromSmarts("[OX2][S](=[O])(=[O])[CH3]")

    def dfs_traverse(node):
        nonlocal alcohol_protection_steps, alcohol_to_leaving_group, leaving_group_substitution

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_str = rsmi.split(">")[0]
                product_str = rsmi.split(">")[-1]

                # Parse reactants and product
                reactants = [Chem.MolFromSmiles(r) for r in reactants_str.split(".") if r]
                product = Chem.MolFromSmiles(product_str) if product_str else None

                if product and any(reactants):
                    # Check for alcohol protection
                    if any(
                        r and r.HasSubstructMatch(Chem.MolFromSmarts("[OH]")) for r in reactants
                    ) and (
                        product.HasSubstructMatch(tbdms_pattern)
                        or product.HasSubstructMatch(benzyl_pattern)
                    ):
                        alcohol_protection_steps += 1
                        print(f"Detected alcohol protection step #{alcohol_protection_steps}")

                    # Check for alcohol to mesylate (leaving group)
                    if any(
                        r and r.HasSubstructMatch(Chem.MolFromSmarts("[OH]")) for r in reactants
                    ) and product.HasSubstructMatch(mesylate_pattern):
                        alcohol_to_leaving_group = True
                        print("Detected alcohol to mesylate conversion")

                    # Check for mesylate to iodide substitution
                    if any(
                        r and r.HasSubstructMatch(mesylate_pattern) for r in reactants
                    ) and product.HasSubstructMatch(Chem.MolFromSmarts("[I]")):
                        leaving_group_substitution = True
                        print("Detected mesylate to iodide substitution")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Check if the strategy is present
    strategy_present = (
        alcohol_protection_steps >= 1 and alcohol_to_leaving_group and leaving_group_substitution
    )

    if strategy_present:
        print("Selective alcohol functionalization strategy detected")

    return strategy_present
