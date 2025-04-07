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
    This function detects a linear synthesis strategy where each step
    has only one main complex reactant (as opposed to convergent synthesis).
    """
    is_linear = True

    def dfs_traverse(node, depth=0):
        nonlocal is_linear

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # Count complex reactants (more than 15 atoms)
                complex_reactants = 0
                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.GetNumAtoms() > 15:  # Arbitrary threshold for "complex"
                        complex_reactants += 1

                if complex_reactants > 1:
                    print(
                        f"Found convergent step at depth {depth} with {complex_reactants} complex reactants"
                    )
                    is_linear = False

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    return is_linear
