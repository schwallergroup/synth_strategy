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
    Detects if the synthesis maintains stereochemistry throughout the route
    by tracking stereochemical notations.
    """
    has_stereocenter = False
    preserves_stereocenter = True

    def dfs_traverse(node):
        nonlocal has_stereocenter, preserves_stereocenter

        if node["type"] == "mol" and "smiles" in node:
            # Check for stereochemical notation in SMILES
            if "@" in node["smiles"]:
                has_stereocenter = True
                print(f"Detected stereocenter in molecule: {node['smiles']}")

        elif (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0]
            products = rsmi.split(">")[-1]

            # If reactants have stereochemistry but products don't, stereochemistry is lost
            if "@" in reactants and "@" not in products:
                preserves_stereocenter = False
                print(
                    f"Stereochemistry lost in reaction at depth {node['metadata'].get('depth', 'unknown')}"
                )

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # Only return True if there is a stereocenter and it's preserved
    result = has_stereocenter and preserves_stereocenter
    print(
        f"Route {'has and preserves stereocenters' if result else 'does not have or does not preserve stereocenters'}"
    )
    return result
