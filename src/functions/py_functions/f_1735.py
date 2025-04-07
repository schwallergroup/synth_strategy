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
    Detects if the synthesis follows a linear pattern where each reaction
    has exactly one product that becomes the reactant for the next step.

    In retrosynthetic analysis, this means each reaction should have exactly
    one non-in_stock reactant molecule.
    """
    is_linear = True

    # Track molecules that are part of the main synthetic pathway
    intermediate_molecules = set()

    def dfs_traverse(node, depth=0):
        nonlocal is_linear

        if node["type"] == "reaction":
            # Get all non-stock reactants (children of type "mol" that are not in_stock)
            non_stock_reactants = [
                child
                for child in node.get("children", [])
                if child["type"] == "mol" and not child.get("in_stock", False)
            ]

            # In a linear synthesis, each reaction in the main pathway should have exactly one non-stock reactant
            # Except for leaf reactions which might have only in-stock reactants
            if len(non_stock_reactants) > 1:
                is_linear = False
                print(
                    f"Non-linear step detected: reaction has {len(non_stock_reactants)} non-stock reactants"
                )

            # Track intermediate molecules for the main pathway
            for child in non_stock_reactants:
                intermediate_molecules.add(child["smiles"])

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # First pass: identify all intermediate molecules
    dfs_traverse(route)

    # Second pass: check if any molecule is used in multiple reactions
    # In a linear synthesis, each intermediate should be used exactly once
    molecule_usage_count = {}

    def count_molecule_usage(node):
        if node["type"] == "mol" and not node.get("in_stock", False):
            smiles = node["smiles"]
            if smiles in molecule_usage_count:
                molecule_usage_count[smiles] += 1
            else:
                molecule_usage_count[smiles] = 1

        # Traverse children
        for child in node.get("children", []):
            count_molecule_usage(child)

    count_molecule_usage(route)

    # Check if any intermediate molecule is used more than once
    # This would indicate a branching or convergent synthesis
    for smiles, count in molecule_usage_count.items():
        if count > 1 and smiles in intermediate_molecules:
            is_linear = False
            print(
                f"Non-linear pattern detected: molecule {smiles} is used in {count} reactions"
            )

    if is_linear:
        print("Linear synthesis pattern detected")

    return is_linear
