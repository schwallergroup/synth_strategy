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
    Detects if the synthesis follows a linear strategy (no convergent steps).
    Linear synthesis is defined as having no more than 2 non-reagent reactants per step.
    """
    is_linear = True

    def dfs_traverse(node):
        nonlocal is_linear

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")

            # Filter out common reagents and solvents
            reagent_patterns = [
                r"O=\[Mn\]",
                r"OH2",
                r"O",
                r"Cl",  # Common reagents
                r"CC\(=O\)O",
                r"CCO",
                r"CCCC",  # Common solvents
            ]

            non_reagent_reactants = []
            for reactant in reactants:
                is_reagent = False
                for pattern in reagent_patterns:
                    if re.search(pattern, reactant):
                        is_reagent = True
                        break
                if not is_reagent:
                    non_reagent_reactants.append(reactant)

            # If more than 2 non-reagent reactants, it's likely convergent
            if len(non_reagent_reactants) > 2:
                is_linear = False
                print(
                    f"Convergent step detected with {len(non_reagent_reactants)} non-reagent reactants"
                )

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return is_linear
