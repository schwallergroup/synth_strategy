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
    Detects if the synthesis follows a linear strategy (each step has only one non-starting material reactant).

    In a linear synthesis, each reaction step should have at most one reactant that is not a starting material.
    Multiple non-starting material reactants indicate a convergent synthesis.

    Args:
        route: A dictionary representing the synthesis route

    Returns:
        bool: True if the synthesis is linear, False otherwise
    """
    is_linear = True

    def dfs_traverse(node):
        nonlocal is_linear

        if node["type"] == "reaction":
            # Get all molecule children (reactants in retrosynthetic direction)
            mol_children = [
                child for child in node.get("children", []) if child["type"] == "mol"
            ]

            # Count non-starting material reactants
            non_starting_reactants = [
                child for child in mol_children if not child.get("in_stock", False)
            ]

            # If more than one non-starting material reactant, it's a convergent synthesis
            if len(non_starting_reactants) > 1:
                is_linear = False

                # Get reaction SMILES for debugging
                reaction_smiles = node.get("metadata", {}).get(
                    "rsmi", "No SMILES available"
                )
                print(
                    f"Found non-linear step (convergent synthesis): {reaction_smiles}"
                )
                print(
                    f"Number of non-starting material reactants: {len(non_starting_reactants)}"
                )

        # Continue traversing the tree
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return is_linear
