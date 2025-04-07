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
    Detects a strategy involving multiple protection/deprotection steps,
    specifically looking for TBDMS protection and benzyl deprotection.
    """
    # Initialize tracking variables
    has_tbdms_protection = False
    has_benzyl_deprotection = False
    protection_deprotection_count = 0

    def dfs_traverse(node):
        nonlocal has_tbdms_protection, has_benzyl_deprotection, protection_deprotection_count

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_str = rsmi.split(">")[0]
                product_str = rsmi.split(">")[-1]

                # Parse reactants and product
                reactants = [Chem.MolFromSmiles(r) for r in reactants_str.split(".") if r]
                product = Chem.MolFromSmiles(product_str) if product_str else None

                if product and any(reactants):
                    # Check for TBDMS protection
                    if any(
                        r
                        and r.HasSubstructMatch(
                            Chem.MolFromSmarts("[Cl][Si]([C])([C])[C]([C])([C])[C]")
                        )
                        for r in reactants
                    ) and product.HasSubstructMatch(
                        Chem.MolFromSmarts("[OX2][Si]([C])([C])[C]([C])([C])[C]")
                    ):
                        has_tbdms_protection = True
                        protection_deprotection_count += 1
                        print("Detected TBDMS protection")

                    # Check for benzyl deprotection
                    if product.HasSubstructMatch(Chem.MolFromSmarts("[OH]")) and any(
                        r
                        and r.HasSubstructMatch(
                            Chem.MolFromSmarts("[OX2][CH2][cX3]1[cX3][cX3][cX3][cX3][cX3]1")
                        )
                        for r in reactants
                    ):
                        has_benzyl_deprotection = True
                        protection_deprotection_count += 1
                        print("Detected benzyl deprotection")

                    # Check for other protection/deprotection steps
                    if (
                        any(
                            r and r.HasSubstructMatch(Chem.MolFromSmarts("[OH]")) for r in reactants
                        )
                        and product.HasSubstructMatch(
                            Chem.MolFromSmarts("[OX2][S](=[O])(=[O])[CH3]")
                        )
                    ) or (
                        any(
                            r
                            and r.HasSubstructMatch(Chem.MolFromSmarts("[OX2][S](=[O])(=[O])[CH3]"))
                            for r in reactants
                        )
                        and product.HasSubstructMatch(Chem.MolFromSmarts("[I]"))
                    ):
                        protection_deprotection_count += 1
                        print("Detected other protection/transformation step")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Check if the strategy is present
    strategy_present = (
        has_tbdms_protection and has_benzyl_deprotection and protection_deprotection_count >= 3
    )

    if strategy_present:
        print(
            f"Protection-deprotection sequence strategy detected with {protection_deprotection_count} steps"
        )

    return strategy_present
