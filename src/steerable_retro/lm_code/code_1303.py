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
    This function detects if the route is a linear synthesis (each reaction has only one product
    that feeds into the next reaction).

    In a linear synthesis:
    1. The final target molecule should have exactly one reaction leading to it
    2. Each intermediate molecule should be the product of exactly one reaction and the reactant for exactly one reaction
    3. Starting materials are only reactants, never products
    4. There should be at least 2 reactions
    """
    # Track all molecules and their roles in reactions
    mol_to_reactions = {}  # Maps molecule SMILES to lists of [producing_rxn, consuming_rxn]
    reaction_count = 0

    def dfs_traverse(node, depth=0):
        nonlocal reaction_count

        if node["type"] == "mol":
            # Initialize tracking for this molecule if not seen before
            smiles = node["smiles"]
            if smiles not in mol_to_reactions:
                mol_to_reactions[smiles] = [None, None]  # [producing_rxn, consuming_rxn]

            # Process children (in retrosynthesis, a mol's children are reactions that produce it)
            for child in node.get("children", []):
                if child["type"] == "reaction":
                    # This reaction produces the current molecule
                    if mol_to_reactions[smiles][0] is not None:
                        print(f"Molecule {smiles} is produced by multiple reactions - not linear")
                        return False
                    mol_to_reactions[smiles][0] = child

                    # Continue traversal
                    if not dfs_traverse(child, depth + 1):
                        return False

        elif node["type"] == "reaction":
            reaction_count += 1

            # Get the product of this reaction (parent molecule)
            product = None
            for parent_node in route.get("children", []):
                if parent_node["type"] == "mol":
                    for rxn_node in parent_node.get("children", []):
                        if rxn_node == node:  # Found the parent molecule
                            product = parent_node["smiles"]
                            break
                    if product:
                        break

            # Process reactants (children of reaction node)
            reactant_count = 0
            for child in node.get("children", []):
                if child["type"] == "mol":
                    reactant_count += 1
                    reactant_smiles = child["smiles"]

                    # Mark this molecule as consumed by this reaction
                    if reactant_smiles not in mol_to_reactions:
                        mol_to_reactions[reactant_smiles] = [None, None]

                    if mol_to_reactions[reactant_smiles][1] is not None:
                        print(
                            f"Molecule {reactant_smiles} is consumed by multiple reactions - not linear"
                        )
                        return False

                    mol_to_reactions[reactant_smiles][1] = node

                    # Continue traversal
                    if not dfs_traverse(child, depth + 1):
                        return False

            # In a linear synthesis, each reaction should have exactly one non-starting material reactant
            # (except possibly the first reaction)
            non_starting_reactants = 0
            for child in node.get("children", []):
                if child["type"] == "mol" and not child.get("in_stock", False):
                    non_starting_reactants += 1

            if non_starting_reactants > 1:
                print(
                    f"Reaction has {non_starting_reactants} non-starting material reactants - not linear"
                )
                return False

        return True

    # Start traversal from the root
    if not dfs_traverse(route):
        return False

    # A linear synthesis should have at least 2 reactions
    if reaction_count < 2:
        print(f"Route has only {reaction_count} reactions - not linear")
        return False

    # Check that we have a single linear path
    # Count molecules that are both produced and consumed
    intermediates = 0
    for smiles, rxns in mol_to_reactions.items():
        if rxns[0] is not None and rxns[1] is not None:
            intermediates += 1

    # In a linear synthesis with N reactions, we should have N-1 intermediates
    if intermediates != reaction_count - 1:
        print(f"Found {intermediates} intermediates but {reaction_count} reactions - not linear")
        return False

    print(f"Route is a linear synthesis with {reaction_count} reactions")
    return True
