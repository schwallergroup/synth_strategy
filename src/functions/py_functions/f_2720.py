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
    Detects if the synthesis route is convergent (multiple branches converge).
    """

    def count_branches(node):
        """Count the number of branches in the synthesis route."""
        if node["type"] == "reaction":
            # Count the number of reactants in this reaction
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                if len(reactants) >= 2:
                    print(
                        f"Convergent synthesis detected with {len(reactants)} reactants"
                    )
                    return True
            except Exception as e:
                print(f"Error analyzing convergent synthesis: {e}")

        # Check children
        for child in node.get("children", []):
            if count_branches(child):
                return True

        return False

    return count_branches(route)
