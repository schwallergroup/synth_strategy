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
    Detects if the synthesis follows a linear strategy without convergent steps.
    """
    is_linear = True

    def count_reactants(node):
        nonlocal is_linear

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # If more than 2 reactants, it might be convergent
                # (allowing for 2 because many reactions have reagents)
                if len(reactants) > 2:
                    is_linear = False
                    print("Found potential convergent step with", len(reactants), "reactants")

        # Traverse children
        for child in node.get("children", []):
            count_reactants(child)

    # Start traversal
    count_reactants(route)

    return is_linear
