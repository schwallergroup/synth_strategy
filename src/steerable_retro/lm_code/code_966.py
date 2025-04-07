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
    This function detects if the synthesis follows a linear fragment assembly pattern.
    """
    fragment_coupling_count = 0
    linear_assembly = True

    def dfs_traverse(node, depth=0):
        nonlocal fragment_coupling_count, linear_assembly

        if node["type"] == "reaction":
            # Get reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]

            # Count number of reactants
            reactant_count = len(reactants_part.split("."))

            if reactant_count > 1:
                fragment_coupling_count += 1

                # If we have more than one fragment coupling at the same depth level,
                # it's not a strictly linear assembly
                if (
                    fragment_coupling_count > 1
                    and node.get("children", [])
                    and any(
                        child["type"] == "reaction"
                        and len(child["metadata"]["rsmi"].split(">")[0].split(".")) > 1
                        for child in node.get("children", [])
                    )
                ):
                    linear_assembly = False

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    return linear_assembly and fragment_coupling_count > 0
