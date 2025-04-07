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
    Detects if the synthesis uses a convergent approach where 3 or more fragments
    are combined in a late-stage step.
    """
    convergent_coupling_detected = False
    late_stage_depth = 1  # Consider depth 0 or 1 as late stage
    min_fragments = 3

    def dfs_traverse(node, depth=0):
        nonlocal convergent_coupling_detected

        if node["type"] == "reaction" and depth <= late_stage_depth:
            # Extract reactants
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")

            # Count number of distinct reactant fragments
            if len(reactants_smiles) >= min_fragments:
                convergent_coupling_detected = True
                print(
                    f"Multi-fragment convergent coupling detected at depth {depth} with {len(reactants_smiles)} fragments"
                )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)
    return convergent_coupling_detected
