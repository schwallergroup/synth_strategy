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
    This function detects if the synthetic route follows a linear (non-convergent) strategy.
    """
    is_linear = True

    def dfs_traverse(node):
        nonlocal is_linear

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # If there are more than 2 reactants, it might be convergent
                # (allowing for 2 because one might be a reagent)
                if len(reactants) > 2:
                    significant_reactants = 0
                    for reactant in reactants:
                        try:
                            mol = Chem.MolFromSmiles(reactant)
                            if (
                                mol and mol.GetNumHeavyAtoms() > 5
                            ):  # Significant molecule, not just a reagent
                                significant_reactants += 1
                        except:
                            continue

                    if significant_reactants > 1:
                        print(
                            "Convergent step detected with multiple significant reactants"
                        )
                        is_linear = False

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return is_linear
