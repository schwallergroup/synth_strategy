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
    This function detects if the synthesis involves a convergent approach where two fragments
    are combined in a late-stage reaction.
    """
    convergent_found = False

    def dfs_traverse(node, depth=0):
        nonlocal convergent_found

        if node["type"] == "reaction" and depth <= 4:  # Focus on late-stage reactions (low depth)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if multiple reactants combine to form a single product
                if len(reactants) >= 2:
                    # Check if reactants are sufficiently complex (not just simple reagents)
                    complex_reactants = 0
                    for reactant in reactants:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and mol.GetNumAtoms() > 6:  # Arbitrary threshold for "complex"
                            complex_reactants += 1

                    if complex_reactants >= 2:
                        print(f"Detected convergent synthesis at depth {depth}")
                        convergent_found = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return convergent_found
