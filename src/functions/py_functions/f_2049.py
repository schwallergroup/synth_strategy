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
    Detects if the synthesis route involves a protection-deprotection sequence,
    particularly focusing on alcohol/phenol protection.
    """
    # Track protection and deprotection events
    protection_events = []
    deprotection_events = []

    def dfs_traverse(node, depth=0):
        nonlocal protection_events, deprotection_events

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Protection patterns
                alcohol_pattern = Chem.MolFromSmarts("[OH]")
                protected_alcohol_patterns = [
                    Chem.MolFromSmarts("[O][CH2][c]"),  # Benzyl ether
                    Chem.MolFromSmarts("[O][Si]"),  # Silyl ether
                    Chem.MolFromSmarts("[O]C(=O)"),  # Ester
                ]

                # Check for protection: alcohol → protected form
                has_alcohol = any(
                    Chem.MolFromSmiles(r) is not None
                    and Chem.MolFromSmiles(r).HasSubstructMatch(alcohol_pattern)
                    for r in reactants
                    if r
                )

                prod_mol = Chem.MolFromSmiles(product) if product else None
                is_protected = prod_mol is not None and any(
                    prod_mol.HasSubstructMatch(pattern)
                    for pattern in protected_alcohol_patterns
                )

                if has_alcohol and is_protected:
                    protection_events.append(depth)
                    print(f"Found protection event at depth {depth}")

                # Check for deprotection: protected form → alcohol
                has_protected = any(
                    Chem.MolFromSmiles(r) is not None
                    and any(
                        Chem.MolFromSmiles(r).HasSubstructMatch(pattern)
                        for pattern in protected_alcohol_patterns
                    )
                    for r in reactants
                    if r
                )

                has_alcohol_product = (
                    prod_mol is not None and prod_mol.HasSubstructMatch(alcohol_pattern)
                )

                if has_protected and has_alcohol_product:
                    deprotection_events.append(depth)
                    print(f"Found deprotection event at depth {depth}")

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if we have both protection and deprotection events
    return len(protection_events) > 0 and len(deprotection_events) > 0
