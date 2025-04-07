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
    This function detects if the synthetic route involves multiple sequential aromatic substitutions.
    """
    aromatic_substitutions = []

    def dfs_traverse(node):
        nonlocal aromatic_substitutions

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                depth = node.get("depth", 0)

                # Check for aromatic chlorination
                if re.search(r"\[Cl:[0-9]+\]", rsmi):
                    aromatic_substitutions.append(("chlorination", depth))
                    print(f"Found aromatic chlorination at depth {depth}: {rsmi}")

                # Check for aromatic nitration
                if re.search(r"\[N\+:[0-9]+\]\(=\[O:[0-9]+\]\)\[O-:[0-9]+\]", rsmi):
                    aromatic_substitutions.append(("nitration", depth))
                    print(f"Found aromatic nitration at depth {depth}: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if we have at least 2 different aromatic substitutions
    substitution_types = set([sub[0] for sub in aromatic_substitutions])
    return len(substitution_types) >= 2
