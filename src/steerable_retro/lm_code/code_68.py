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
    Detects if the synthesis maintains the core scaffold while modifying functional groups.
    This is characterized by preserving aromatic systems throughout the synthesis.
    """
    aromatic_systems = []
    scaffold_preserved = True

    def dfs_traverse(node):
        nonlocal aromatic_systems, scaffold_preserved

        if node["type"] == "mol":
            smiles = node["smiles"]
            if smiles:
                # Count aromatic rings
                aromatic_count = smiles.count("c")
                if aromatic_count > 0:
                    if not aromatic_systems:
                        aromatic_systems.append(aromatic_count)
                    else:
                        # If the number of aromatic carbons changes significantly,
                        # the scaffold is not preserved
                        if abs(aromatic_count - aromatic_systems[0]) > 2:
                            print(
                                f"Scaffold changed: aromatic count from {aromatic_systems[0]} to {aromatic_count}"
                            )
                            scaffold_preserved = False

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return scaffold_preserved and len(aromatic_systems) > 0
