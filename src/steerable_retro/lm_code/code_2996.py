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
    Detects if the synthesis follows a convergent strategy with at least
    two complex fragments combined in late-stage reactions (depth 0 or 1).
    """
    result = False

    def dfs_traverse(node, depth=0):
        nonlocal result

        if node["type"] == "reaction" and depth <= 1:
            # Extract reactants
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")

            # Check if we have at least 2 complex reactants
            complex_reactants = 0
            for reactant in reactants_smiles:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol:
                        # Define complexity as having at least 10 heavy atoms
                        if mol.GetNumHeavyAtoms() >= 10:
                            complex_reactants += 1
                except:
                    continue

            if complex_reactants >= 2:
                print(
                    f"Detected convergent synthesis at depth {depth} with {complex_reactants} complex fragments"
                )
                result = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return result
