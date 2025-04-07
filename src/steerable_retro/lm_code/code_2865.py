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
    This function detects if the synthesis follows a linear strategy (no convergent steps).
    """
    is_linear = True

    def dfs_traverse(node, depth=0):
        nonlocal is_linear

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # Count non-reagent reactants (simplified approach)
                significant_reactants = 0
                for reactant in reactants:
                    # Skip common reagents (simplified)
                    if (
                        "O=P(Cl)(Cl)Cl" in reactant
                        or "O=[N+]([O-])[O-]" in reactant
                        or "N" == reactant
                    ):
                        continue
                    significant_reactants += 1

                if significant_reactants > 2:
                    print(
                        f"Found convergent step with {significant_reactants} significant reactants at depth {depth}"
                    )
                    is_linear = False

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return is_linear
