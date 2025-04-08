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
    Detects if the synthesis preserves stereochemistry throughout the route
    """
    has_stereocenter_in_early_stage = False
    maintains_stereocenter_to_final = False

    def dfs_traverse(node, depth=0, current_path=None):
        nonlocal has_stereocenter_in_early_stage, maintains_stereocenter_to_final

        if current_path is None:
            current_path = []

        if node["type"] == "mol":
            # Check for stereocenters
            if "smiles" in node:
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol:
                    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=False)
                    if len(chiral_centers) > 0:
                        if depth > 2:  # Early stage
                            has_stereocenter_in_early_stage = True
                            print(f"Found stereocenter in early stage (depth {depth})")
                        if depth == 0:  # Final product
                            maintains_stereocenter_to_final = True
                            print("Final product has stereocenter")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1, current_path + [node])

    # Start traversal from root
    dfs_traverse(route)
    return has_stereocenter_in_early_stage and maintains_stereocenter_to_final
