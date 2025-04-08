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
    Detects a convergent synthesis strategy with multiple fragment couplings.
    Counts reactions with 2 or more reactants.
    """
    fragment_coupling_count = 0

    def dfs_traverse(node):
        nonlocal fragment_coupling_count

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # Check if we have multiple reactants (convergent)
                if len(reactants) >= 2:
                    # Check if reactants are complex enough (not just simple reagents)
                    complex_reactants = 0
                    for reactant in reactants:
                        try:
                            mol = Chem.MolFromSmiles(reactant)
                            if mol and mol.GetNumAtoms() > 6:  # Arbitrary threshold for "complex"
                                complex_reactants += 1
                        except:
                            continue

                    if complex_reactants >= 2:
                        fragment_coupling_count += 1
                        print(
                            f"Found fragment coupling reaction with {complex_reactants} complex reactants"
                        )

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return fragment_coupling_count >= 2  # At least 2 fragment couplings for convergent strategy
