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
    This function detects a protection-deprotection strategy where a phenol is protected
    as an allyl ether and later deprotected back to phenol.
    """
    # Track if we find both protection and deprotection
    protection_found = False
    deprotection_found = False

    # SMARTS patterns
    phenol_pattern = Chem.MolFromSmarts("[OH]-c")
    allyl_ether_pattern = Chem.MolFromSmarts("[O]-[CH2]-[CH]=[CH2]")

    def dfs_traverse(node):
        nonlocal protection_found, deprotection_found

        if node["type"] == "reaction":
            # Extract reactants and products
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            products_part = rsmi.split(">")[-1]

            reactants = [Chem.MolFromSmiles(r) for r in reactants_part.split(".") if r]
            product = Chem.MolFromSmiles(products_part)

            # Check for protection (phenol → allyl ether)
            if (
                any(
                    r is not None and r.HasSubstructMatch(phenol_pattern)
                    for r in reactants
                )
                and product is not None
                and product.HasSubstructMatch(allyl_ether_pattern)
            ):
                protection_found = True
                print("Found phenol protection with allyl group")

            # Check for deprotection (allyl ether → phenol)
            if (
                any(
                    r is not None and r.HasSubstructMatch(allyl_ether_pattern)
                    for r in reactants
                )
                and product is not None
                and product.HasSubstructMatch(phenol_pattern)
            ):
                deprotection_found = True
                print("Found allyl ether deprotection to phenol")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if both protection and deprotection were found
    return protection_found and deprotection_found
