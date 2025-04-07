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
    Detects if the synthesis follows a linear strategy (as opposed to convergent).
    Checks if most reactions have only one non-commercial reactant.
    """
    reaction_count = 0
    linear_reaction_count = 0

    def dfs_traverse(node):
        nonlocal reaction_count, linear_reaction_count

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                reaction_count += 1

                # Count non-commercial reactants
                non_commercial_count = 0
                for child in node.get("children", []):
                    if child["type"] == "mol" and not child.get("in_stock", False):
                        non_commercial_count += 1

                if non_commercial_count <= 1:
                    linear_reaction_count += 1
                    print(
                        f"Found linear reaction step with {non_commercial_count} non-commercial reactants"
                    )

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # If at least 75% of reactions are linear, consider it a linear synthesis
    if reaction_count > 0:
        linear_ratio = linear_reaction_count / reaction_count
        print(f"Linear reaction ratio: {linear_ratio} ({linear_reaction_count}/{reaction_count})")
        return linear_ratio >= 0.75

    return False
