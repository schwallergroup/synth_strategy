#!/bin/python

"""LM-defined function for strategy description."""

from rdkit.Chem import AllChem, rdFMCS
import copy
from collections import deque
import rdkit.Chem as Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdChemReactions
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
import rdkit.Chem.rdFMCS
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem import AllChem, Descriptors, Lipinski
from rdkit.Chem import rdmolops
import re
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import AllChem, Descriptors
import traceback
import rdkit
from collections import Counter


def main(route):
    """
    Detects if the synthesis follows a linear strategy (each reaction has only one
    non-starting material reactant) rather than a convergent strategy.

    A linear synthesis builds complexity by sequentially adding starting materials to a growing molecule.
    A convergent synthesis combines complex intermediates in later stages.

    Args:
        route: A dictionary representing the synthesis route following the SynthesisRoute schema

    Returns:
        bool: True if the synthesis is linear, False if it's convergent
    """
    is_linear = True

    def dfs_traverse(node, depth=0):
        nonlocal is_linear

        # Process reaction nodes to check for convergent steps
        if node["type"] == "reaction":
            # In retrosynthetic traversal, we're looking at reactants that form the product
            # For a linear synthesis, at most one of these reactants should be a complex intermediate
            complex_intermediates = []

            for child in node.get("children", []):
                if child["type"] == "mol" and not child.get("in_stock", False):
                    complex_intermediates.append(child)

            # If more than one complex intermediate is used in this reaction, it's convergent
            if len(complex_intermediates) > 1:
                # Check if these are truly separate synthetic pathways
                # Sometimes multiple "non-starting materials" are actually part of the same pathway
                # and don't represent a convergent synthesis

                # For this implementation, we'll consider it convergent if there are multiple
                # complex intermediates that each have their own reaction nodes as children
                separate_pathways = 0
                for intermediate in complex_intermediates:
                    has_reaction_children = False
                    for child in intermediate.get("children", []):
                        if child["type"] == "reaction":
                            has_reaction_children = True
                            break
                    if has_reaction_children:
                        separate_pathways += 1

                if separate_pathways > 1:
                    print(
                        f"Depth {depth}: Found convergent step with {separate_pathways} separate synthetic pathways"
                    )
                    intermediate_smiles = [
                        inter.get("smiles", "unknown") for inter in complex_intermediates
                    ]
                    print(f"Complex intermediates: {intermediate_smiles}")
                    is_linear = False

        # Continue traversing the synthesis tree
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root node
    dfs_traverse(route)

    print(f"Synthesis strategy: {'Linear' if is_linear else 'Convergent'}")
    return is_linear
