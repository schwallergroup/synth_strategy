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
    This function detects if the synthesis route preserves stereochemistry throughout.
    """
    stereocenters_by_depth = {}

    def dfs_traverse(node, depth=0):
        if node["type"] == "mol":
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol is not None:
                # Count chiral centers with defined stereochemistry
                chiral_centers = 0
                for atom in mol.GetAtoms():
                    if atom.HasProp("_CIPCode"):  # Atom has defined stereochemistry
                        chiral_centers += 1

                stereocenters_by_depth[depth] = chiral_centers
                if chiral_centers > 0:
                    print(f"Found {chiral_centers} stereocenters at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if stereocenters are present throughout the synthesis
    depths_with_stereo = [
        depth for depth, count in stereocenters_by_depth.items() if count > 0
    ]

    # If we have stereocenters at multiple depths including the final product (depth 0)
    preserves_stereochemistry = 0 in depths_with_stereo and len(depths_with_stereo) > 1
    print(f"Stereocenter preservation: {preserves_stereochemistry}")
    return preserves_stereochemistry
