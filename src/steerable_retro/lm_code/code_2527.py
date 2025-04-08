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
    Detects if the synthesis follows a linear functionalization strategy
    rather than a convergent approach.
    """
    # Track the maximum branching factor in the synthesis tree
    max_branching = 0

    def dfs_traverse(node):
        nonlocal max_branching

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # Count significant reactants (more than just reagents)
                significant_reactants = 0
                for reactant in reactants:
                    r_mol = Chem.MolFromSmiles(reactant)
                    if (
                        r_mol and r_mol.GetNumAtoms() > 3
                    ):  # Arbitrary threshold to exclude small reagents
                        significant_reactants += 1

                max_branching = max(max_branching, significant_reactants)

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # If max_branching is consistently low (≤2), it suggests a linear strategy
    is_linear = max_branching <= 2
    print(f"Maximum branching factor: {max_branching}, Linear strategy: {is_linear}")
    return is_linear
