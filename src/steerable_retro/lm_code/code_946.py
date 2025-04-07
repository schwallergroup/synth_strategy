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
    Detects a strategy involving fragment assembly through heteroatom alkylations (N and O).
    Looks for both N-alkylation and O-alkylation in the synthesis route.
    """
    # Track if we've found both types of alkylations
    n_alkylation_found = False
    o_alkylation_found = False

    def dfs_traverse(node, depth=0):
        nonlocal n_alkylation_found, o_alkylation_found

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for N-alkylation: N-H + X-C → N-C
            # Look for secondary amine in reactants and tertiary amine in product
            if (
                any(
                    Chem.MolFromSmiles(r)
                    and Chem.MolFromSmiles(r).HasSubstructMatch(
                        Chem.MolFromSmarts("[N;H1]([C])[C]")
                    )
                    for r in reactants
                )
                and any(
                    Chem.MolFromSmiles(r)
                    and Chem.MolFromSmiles(r).HasSubstructMatch(Chem.MolFromSmarts("[C][Cl,Br,I]"))
                    for r in reactants
                )
                and Chem.MolFromSmiles(product)
                and Chem.MolFromSmiles(product).HasSubstructMatch(
                    Chem.MolFromSmarts("[N]([C])([C])[C]")
                )
            ):
                n_alkylation_found = True
                print(f"Found N-alkylation at depth {depth}")

            # Check for O-alkylation: O-H + X-C → O-C
            if (
                any(
                    Chem.MolFromSmiles(r)
                    and Chem.MolFromSmiles(r).HasSubstructMatch(Chem.MolFromSmarts("[O;H1][C]"))
                    for r in reactants
                )
                and any(
                    Chem.MolFromSmiles(r)
                    and Chem.MolFromSmiles(r).HasSubstructMatch(Chem.MolFromSmarts("[C][Cl,Br,I]"))
                    for r in reactants
                )
                and Chem.MolFromSmiles(product)
                and Chem.MolFromSmiles(product).HasSubstructMatch(Chem.MolFromSmarts("[O]([C])[C]"))
            ):
                o_alkylation_found = True
                print(f"Found O-alkylation at depth {depth}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Return True if both types of alkylations are found
    return n_alkylation_found and o_alkylation_found
