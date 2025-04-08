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
    Detects a linear fragment assembly strategy where the molecule is built
    sequentially without convergent steps.
    """
    # Track the maximum number of reactants in any reaction
    max_reactants = 0
    has_coupling_reaction = False

    def dfs_traverse(node):
        nonlocal max_reactants, has_coupling_reaction

        if node["type"] == "reaction":
            # Extract reactants
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            reactants = reactants_part.split(".")

            # Count number of significant reactants (excluding small molecules/reagents)
            significant_reactants = 0
            for r in reactants:
                # Skip empty strings or very small molecules (likely reagents)
                if r and len(r) > 5:  # Arbitrary threshold to exclude small reagents
                    significant_reactants += 1

            max_reactants = max(max_reactants, significant_reactants)

            # Check if this is a coupling reaction (joining two significant fragments)
            if significant_reactants >= 2:
                has_coupling_reaction = True
                print(
                    f"Detected coupling reaction with {significant_reactants} significant reactants"
                )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # A linear strategy typically has at most 2 significant reactants in any step
    # and should have at least one coupling reaction
    is_linear = max_reactants <= 2 and has_coupling_reaction

    if is_linear:
        print("Linear fragment assembly strategy detected")
        return True
    else:
        print("Linear fragment assembly strategy not detected")
        return False
