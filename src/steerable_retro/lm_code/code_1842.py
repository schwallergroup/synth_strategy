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
    Checks if most reactions have only 1-2 reactants.
    """
    reaction_count = 0
    convergent_count = 0

    def dfs_traverse(node):
        nonlocal reaction_count, convergent_count

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")

            reaction_count += 1
            # Count complex reactants (not simple reagents)
            complex_reactants = 0
            for reactant in reactants:
                # Simple heuristic: complex molecules have more than 10 atoms
                if len(reactant) > 10 and not re.search(r"[BPOS]|Cl|Br|[0-9]", reactant):
                    complex_reactants += 1

            if complex_reactants > 1:
                convergent_count += 1
                print(f"Convergent step detected: {rsmi}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # If less than 20% of reactions are convergent, consider it linear
    return reaction_count > 0 and convergent_count / reaction_count < 0.2
