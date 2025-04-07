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
    This function detects if the synthesis follows a convergent approach
    by identifying steps with multiple reactants that aren't reagents.
    """
    convergent_steps = 0

    def dfs_traverse(node):
        nonlocal convergent_steps

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # Count significant reactants (not small reagents)
                significant_reactants = 0
                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if (
                        mol and mol.GetNumHeavyAtoms() > 6
                    ):  # Arbitrary threshold to exclude small reagents
                        significant_reactants += 1

                if significant_reactants >= 2:
                    convergent_steps += 1
                    print(
                        f"Convergent step detected with {significant_reactants} significant reactants"
                    )

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return convergent_steps > 0
