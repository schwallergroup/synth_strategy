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
    This function detects if the synthesis follows a linear strategy
    (each reaction has only one product that serves as a reactant in the next step).

    In a retrosynthetic tree:
    - The root is the target molecule
    - Reaction nodes have molecule children (reactants in forward direction)
    - A linear synthesis has a single path from root to leaves with each reaction having exactly one non-terminal reactant
    """
    # Track if synthesis is linear and count reactions
    is_linear = True
    reaction_count = 0

    # Check if the root is a molecule
    if route["type"] != "mol":
        print("Root node is not a molecule")
        return False

    current_node = route

    # Follow the path from target molecule through the synthesis
    while current_node.get("children", []):
        # If current node is a molecule, it should have exactly one reaction child
        if current_node["type"] == "mol":
            if len(current_node.get("children", [])) != 1:
                print(
                    f"Non-linear step detected: molecule has {len(current_node.get('children', []))} reaction children"
                )
                is_linear = False
                break

            # Move to the reaction child
            current_node = current_node["children"][0]

        # If current node is a reaction
        elif current_node["type"] == "reaction":
            reaction_count += 1

            # Count non-terminal reactants (molecules that aren't starting materials)
            non_terminal_reactants = 0
            next_mol = None

            # Count non-terminal reactants and find the next molecule in the path
            for child in current_node.get("children", []):
                if child["type"] == "mol" and not child.get("in_stock", False):
                    non_terminal_reactants += 1
                    next_mol = child

            # In a linear synthesis, a reaction should have at most one non-terminal reactant
            if non_terminal_reactants > 1:
                print(
                    f"Non-linear step detected: reaction has {non_terminal_reactants} non-terminal reactants"
                )
                is_linear = False
                break
            elif non_terminal_reactants == 0:
                # This is a terminal reaction (all reactants are in_stock)
                print(f"Reached terminal reaction (all reactants are in_stock)")
                break

            # Move to the next molecule in the path
            current_node = next_mol

        else:
            print(f"Unknown node type: {current_node['type']}")
            is_linear = False
            break

    # A linear synthesis should have at least 2 reactions
    if reaction_count < 2:
        print(f"Not enough reactions: {reaction_count} (minimum 2 required)")
        is_linear = False

    if is_linear:
        print(f"Linear synthesis strategy detected with {reaction_count} reactions")

    return is_linear
