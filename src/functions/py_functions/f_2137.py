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
    Detects a convergent synthesis where multiple fragments are prepared separately
    and then coupled together.
    """
    fragment_counts = {}  # depth -> count of fragments
    coupling_reactions = []  # list of (depth, num_reactants)

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")

            # Count number of reactants
            num_reactants = len(reactants)

            # Record reactions with multiple reactants (potential coupling reactions)
            if num_reactants > 1:
                coupling_reactions.append((depth, num_reactants))
                print(
                    f"Found potential coupling reaction at depth {depth} with {num_reactants} reactants"
                )

            # Update fragment count at this depth
            if depth in fragment_counts:
                fragment_counts[depth] += (
                    num_reactants - 1
                )  # -1 because we're counting additional fragments
            else:
                fragment_counts[depth] = num_reactants - 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Analyze the fragment counts and coupling reactions
    total_fragments = sum(fragment_counts.values())
    has_multi_fragment_coupling = any(
        num_reactants >= 2 for _, num_reactants in coupling_reactions
    )

    # Check if this is a convergent synthesis (multiple fragments, coupling reactions)
    is_convergent = total_fragments >= 2 and has_multi_fragment_coupling

    print(f"Fragment counts by depth: {fragment_counts}")
    print(f"Coupling reactions: {coupling_reactions}")
    print(f"Total fragments: {total_fragments}")
    print(f"Convergent fragment coupling strategy detected: {is_convergent}")

    return is_convergent
