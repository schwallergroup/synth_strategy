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
    This function detects a linear synthesis strategy with a late-stage fragment coupling,
    where two significant fragments are combined in the final steps.
    """
    late_fragment_coupling = False

    def dfs_traverse(node, depth=0):
        nonlocal late_fragment_coupling

        if node["type"] == "reaction" and depth <= 1:  # Focus on late-stage reactions
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # Check if we have at least two significant fragments
                if len(reactants) >= 2:
                    significant_fragments = 0

                    for reactant in reactants:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol:
                            # Consider a fragment significant if it has at least 6 atoms
                            if mol.GetNumAtoms() >= 6:
                                significant_fragments += 1

                    if significant_fragments >= 2:
                        late_fragment_coupling = True
                        print("Late-stage fragment coupling detected at depth", depth)

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    return late_fragment_coupling
