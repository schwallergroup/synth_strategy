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
    Detects if the synthesis follows a linear strategy without convergent steps.

    A linear synthesis is one where each reaction node has at most one child that is a reaction node,
    meaning there are no convergent steps where multiple complex intermediates come together.
    """
    is_linear = True

    def dfs_traverse(node, depth=0):
        nonlocal is_linear

        # Skip processing if we already know it's not linear
        if not is_linear:
            return

        if node["type"] == "reaction":
            print(
                f"Checking reaction at depth {depth}: {node.get('metadata', {}).get('rsmi', 'No RSMI')}"
            )

            # Count reaction children (non-terminal nodes)
            reaction_children = 0
            for child in node.get("children", []):
                if child["type"] == "mol" and not child.get("in_stock", False):
                    # This is an intermediate molecule, count its reaction children
                    mol_reaction_children = sum(
                        1
                        for grandchild in child.get("children", [])
                        if grandchild["type"] == "reaction"
                    )
                    reaction_children += mol_reaction_children

            # If more than one reaction path continues from this node, it's convergent
            if reaction_children > 1:
                print(
                    f"Found convergent structure at depth {depth} with {reaction_children} reaction paths"
                )
                is_linear = False

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return is_linear
