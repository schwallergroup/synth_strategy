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
    Detects if the synthesis preserves a stereocenter throughout the route.
    """
    stereocenter_in_final = False
    stereocenter_in_intermediate = False

    def dfs_traverse(node, depth=0):
        nonlocal stereocenter_in_final, stereocenter_in_intermediate

        if node["type"] == "mol" and "smiles" in node:
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Find chiral centers
                chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=False)

                if len(chiral_centers) > 0:
                    if depth == 0:
                        stereocenter_in_final = True
                    else:
                        stereocenter_in_intermediate = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    print(f"Stereocenter in final product: {stereocenter_in_final}")
    print(f"Stereocenter in intermediate: {stereocenter_in_intermediate}")

    return stereocenter_in_final and stereocenter_in_intermediate
