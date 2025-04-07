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
    This function detects if the synthesis follows a linear strategy where each step
    builds on a single product from the previous step (as opposed to convergent synthesis).
    """
    is_linear = True
    reaction_count = 0

    def dfs_traverse(node, depth=0):
        nonlocal is_linear, reaction_count

        if node["type"] == "reaction":
            reaction_count += 1
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # Count non-empty reactants
                reactant_count = sum(
                    1 for r in reactants if r and not (r.strip() in ["Br2", "HBr", "NBS"])
                )

                # If more than 2 substantial reactants, it's likely not linear
                if reactant_count > 2:
                    print(
                        f"Non-linear synthesis detected at depth {depth} with {reactant_count} reactants"
                    )
                    is_linear = False

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    # Must have at least 3 reactions to be considered a meaningful linear synthesis
    return is_linear and reaction_count >= 3
