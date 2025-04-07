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
    Detects if the synthesis follows a convergent approach where two complex fragments
    are prepared separately and then coupled in a late-stage reaction.
    """
    # Track if we've found a convergent step
    convergent_step_found = False

    def dfs_traverse(node, depth=0):
        nonlocal convergent_step_found

        if node["type"] == "reaction" and depth <= 1:  # Focus on late-stage reactions
            # Extract reactants
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")

            # Check if we have at least 2 complex reactants
            complex_reactants = 0
            for reactant in reactants_smiles:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol:
                        # Define complexity as having more than 10 atoms
                        if mol.GetNumAtoms() > 10:
                            complex_reactants += 1
                except:
                    continue

            if complex_reactants >= 2:
                convergent_step_found = True
                print(
                    f"Convergent synthesis detected at depth {depth} with {complex_reactants} complex reactants"
                )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return convergent_step_found
