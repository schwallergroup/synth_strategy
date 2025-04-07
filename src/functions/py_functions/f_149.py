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
    Detects if the synthesis follows a linear strategy (as opposed to convergent).
    """
    max_fragments_per_step = 0

    def dfs_traverse(node):
        nonlocal max_fragments_per_step

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")

            # Count non-trivial reactants (exclude simple reagents)
            complex_reactants = 0
            for reactant in reactants:
                # Simple reagents typically have fewer atoms
                mol = Chem.MolFromSmiles(reactant)
                if mol and mol.GetNumAtoms() > 6:  # Arbitrary threshold for "complex"
                    complex_reactants += 1

            max_fragments_per_step = max(max_fragments_per_step, complex_reactants)

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # Linear synthesis typically has at most 2 complex fragments per step
    is_linear = max_fragments_per_step <= 2
    print(f"Maximum complex fragments per step: {max_fragments_per_step}")
    print(f"Synthesis strategy is {'linear' if is_linear else 'convergent'}")

    return is_linear
