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
    This function detects if the synthesis follows a linear strategy without convergent steps.
    Linear synthesis means each reaction has only one non-starting material reactant.
    """
    is_linear = True

    def dfs_traverse(node):
        nonlocal is_linear

        if node["type"] == "reaction":
            # Check if this is a linear step (only one non-starting material reactant)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                reactants = reactants_part.split(".")

                # Count non-starting material reactants
                non_starting_material_count = 0
                for child in node.get("children", []):
                    if child["type"] == "mol" and not child.get("in_stock", False):
                        non_starting_material_count += 1

                if non_starting_material_count > 1:
                    is_linear = False
                    print(
                        f"Found non-linear step with {non_starting_material_count} non-starting material reactants"
                    )

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(f"Synthesis is {'linear' if is_linear else 'convergent'}")
    return is_linear
