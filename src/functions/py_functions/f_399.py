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
    This function detects if the synthetic route involves protecting group strategy,
    particularly Boc protection/deprotection.
    """
    has_protection_deprotection = False
    protection_events = 0
    deprotection_events = 0

    def dfs_traverse(node):
        nonlocal has_protection_deprotection, protection_events, deprotection_events

        if node["type"] == "reaction":
            # Get reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for Boc group
            boc_pattern = Chem.MolFromSmarts("[CH3]C([CH3])([CH3])[O]C(=[O])[NH]")

            # Check if Boc is added (protection)
            product_mol = Chem.MolFromSmiles(product_smiles)
            has_boc_in_product = (
                product_mol.HasSubstructMatch(boc_pattern) if product_mol else False
            )

            has_boc_in_reactants = False
            for r_smiles in reactants_smiles:
                r_mol = Chem.MolFromSmiles(r_smiles)
                if r_mol and r_mol.HasSubstructMatch(boc_pattern):
                    has_boc_in_reactants = True
                    break

            # Protection: Boc not in reactants but in product
            if not has_boc_in_reactants and has_boc_in_product:
                protection_events += 1
                print("Boc protection detected")

            # Deprotection: Boc in reactants but not in product
            if has_boc_in_reactants and not has_boc_in_product:
                deprotection_events += 1
                print("Boc deprotection detected")

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from root
    dfs_traverse(route)

    # Consider it a protecting group strategy if there's at least one protection or deprotection event
    has_protection_deprotection = protection_events > 0 or deprotection_events > 0
    if has_protection_deprotection:
        print(
            f"Protecting group strategy detected: {protection_events} protection events, {deprotection_events} deprotection events"
        )

    return has_protection_deprotection
