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
    This function detects if the synthesis employs a convergent fragment coupling strategy
    with at least 3 distinct coupling events.
    """
    coupling_events = 0

    def dfs_traverse(node):
        nonlocal coupling_events

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")

            # If there are multiple reactants, this might be a coupling event
            if len(reactants) >= 2:
                # Check if reactants are complex enough (not just simple reagents)
                complex_reactants = 0
                for r in reactants:
                    mol = Chem.MolFromSmiles(r)
                    if mol and mol.GetNumAtoms() > 5:  # Arbitrary threshold for "complex" molecules
                        complex_reactants += 1

                if complex_reactants >= 2:
                    coupling_events += 1
                    print(f"Found coupling event: {rsmi}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return coupling_events >= 3
