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
    This function detects if the synthesis follows a linear strategy (vs convergent).
    Linear synthesis typically has only 1-2 reactants per step.
    """
    is_linear = True

    def dfs_traverse(node):
        nonlocal is_linear

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")

            # Count non-trivial reactants (excluding simple reagents)
            complex_reactants = 0
            for reactant in reactants:
                mol = Chem.MolFromSmiles(reactant)
                if mol and mol.GetNumAtoms() > 5:  # Arbitrary threshold for "complex"
                    complex_reactants += 1

            if complex_reactants > 2:
                is_linear = False
                print(
                    f"Non-linear step detected with {complex_reactants} complex reactants: {rsmi}"
                )

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return is_linear
