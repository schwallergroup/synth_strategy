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


def main(route, max_reactants_per_step=2):
    """
    This function detects if the route follows a linear synthesis strategy with limited reactants per step.
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
            reactants_smiles = rsmi.split(">")[0].split(".")

            # Count non-empty reactants
            reactant_count = sum(1 for r in reactants_smiles if r)

            if reactant_count > max_reactants_per_step:
                is_linear = False
                print(
                    f"Non-linear step detected with {reactant_count} reactants: {rsmi}"
                )

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(f"Linear synthesis strategy: {is_linear}")
    return is_linear
