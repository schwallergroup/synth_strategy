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
    This function detects if the synthetic route follows a linear strategy (no convergent steps).
    """
    is_linear = True

    def dfs_traverse(node):
        nonlocal is_linear

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # If there are more than 2 reactants, it might be a convergent step
                # (allowing for 2 because many reactions have a reagent alongside the main reactant)
                if len(reactants) > 2:
                    # Check if at least 2 reactants are complex (not simple reagents)
                    complex_reactants = 0
                    for reactant in reactants:
                        try:
                            mol = Chem.MolFromSmiles(reactant)
                            if mol and mol.GetNumAtoms() > 5:  # Arbitrary threshold for "complex"
                                complex_reactants += 1
                        except:
                            continue

                    if complex_reactants >= 2:
                        print(
                            "Convergent step detected with", complex_reactants, "complex reactants"
                        )
                        is_linear = False

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return is_linear
