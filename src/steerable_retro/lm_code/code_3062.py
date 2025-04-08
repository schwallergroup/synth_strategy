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
    Detects the use of a stereocentered building block that is carried through
    multiple steps of the synthesis.
    """
    stereocenters_by_depth = {}

    def dfs_traverse(node, depth=0):
        if node["type"] == "mol":
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Find stereocenters
                chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=False)
                if chiral_centers:
                    if depth not in stereocenters_by_depth:
                        stereocenters_by_depth[depth] = []
                    stereocenters_by_depth[depth].append(len(chiral_centers))
                    print(f"Found {len(chiral_centers)} stereocenters at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if stereocenters are present in at least 3 different depths
    if len(stereocenters_by_depth) >= 3:
        return True

    return False
