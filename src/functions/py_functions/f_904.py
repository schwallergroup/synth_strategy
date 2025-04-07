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
    This function detects if the synthesis follows a convergent approach
    by counting the number of fragments combined in the synthesis.
    """
    max_fragments = 0

    def count_fragments(rsmi):
        reactants = rsmi.split(">")[0].split(".")
        return len(reactants)

    def dfs_traverse(node):
        nonlocal max_fragments

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            fragments = count_fragments(node["metadata"]["rsmi"])
            if fragments > max_fragments:
                max_fragments = fragments

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(f"Maximum fragments in a single reaction: {max_fragments}")
    return (
        max_fragments >= 2
    )  # Consider convergent if at least 2 fragments are combined
