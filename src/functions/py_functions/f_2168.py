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
    Detects orthogonal protection strategy where multiple alcohols in the same molecule
    are protected with different groups (TBDMS and benzyl), followed by selective
    deprotection and functional group conversion.
    """
    # Initialize tracking variables
    has_tbdms_protection = False
    has_benzyl_protection = False
    has_benzyl_deprotection = False
    has_alcohol_to_mesylate = False
    has_mesylate_to_iodide = False

    # SMARTS patterns
    tbdms_pattern = Chem.MolFromSmarts("[OX2][Si]([C])([C])[C]([C])([C])[C]")
    benzyl_pattern = Chem.MolFromSmarts("[OX2][CH2][cX3]1[cX3][cX3][cX3][cX3][cX3]1")
    mesylate_pattern = Chem.MolFromSmarts("[OX2][S](=[O])(=[O])[CH3]")

    def dfs_traverse(node):
        nonlocal has_tbdms_protection, has_benzyl_protection, has_benzyl_deprotection
        nonlocal has_alcohol_to_mesylate, has_mesylate_to_iodide

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_str = rsmi.split(">")[0]
                product_str = rsmi.split(">")[-1]

                # Parse reactants and product
                reactants = [
                    Chem.MolFromSmiles(r) for r in reactants_str.split(".") if r
                ]
                product = Chem.MolFromSmiles(product_str) if product_str else None

                if product and any(reactants):
                    # Check for TBDMS protection
                    if any(
                        r
                        and r.HasSubstructMatch(
                            Chem.MolFromSmarts("[Cl][Si]([C])([C])[C]([C])([C])[C]")
                        )
                        for r in reactants
                    ) and product.HasSubstructMatch(tbdms_pattern):
                        has_tbdms_protection = True
                        print("Detected TBDMS protection")

                    # Check for benzyl protection/deprotection
                    if any(
                        r
                        and r.HasSubstructMatch(
                            Chem.MolFromSmarts("[c]1[c][c][c][c][c]1[CH2][O]")
                        )
                        for r in reactants
                    ) and product.HasSubstructMatch(Chem.MolFromSmarts("[OH]")):
                        has_benzyl_deprotection = True
                        print("Detected benzyl deprotection")

                    # Check for alcohol to mesylate conversion
                    if (
                        any(
                            r and r.HasSubstructMatch(Chem.MolFromSmarts("[OH]"))
                            for r in reactants
                        )
                        and any(
                            r
                            and r.HasSubstructMatch(
                                Chem.MolFromSmarts("[Cl][S](=[O])(=[O])[CH3]")
                            )
                            for r in reactants
                        )
                        and product.HasSubstructMatch(mesylate_pattern)
                    ):
                        has_alcohol_to_mesylate = True
                        print("Detected alcohol to mesylate conversion")

                    # Check for mesylate to iodide conversion
                    if product.HasSubstructMatch(Chem.MolFromSmarts("[I]")) and any(
                        r and r.HasSubstructMatch(mesylate_pattern) for r in reactants
                    ):
                        has_mesylate_to_iodide = True
                        print("Detected mesylate to iodide conversion")

                    # Check for benzyl protection
                    benzyl_in_product = (
                        product.HasSubstructMatch(benzyl_pattern) if product else False
                    )
                    if benzyl_in_product:
                        has_benzyl_protection = True
                        print("Detected benzyl protection")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Check if the strategy is present
    strategy_present = (
        has_tbdms_protection
        and (has_benzyl_protection or has_benzyl_deprotection)
        and (has_alcohol_to_mesylate or has_mesylate_to_iodide)
    )

    if strategy_present:
        print("Orthogonal alcohol protection strategy detected")

    return strategy_present
