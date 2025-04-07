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
    This function detects a convergent synthesis strategy with multiple fragment couplings.
    """
    fragment_couplings = 0

    def dfs_traverse(node, depth=0):
        nonlocal fragment_couplings

        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for fragment coupling (multiple reactants -> single product)
                if len(reactants) >= 2 and product:
                    # Count non-trivial reactants (exclude small molecules/reagents)
                    significant_reactants = 0
                    for r in reactants:
                        r_mol = Chem.MolFromSmiles(r)
                        if (
                            r_mol and r_mol.GetNumHeavyAtoms() > 6
                        ):  # Arbitrary threshold for "significant" fragment
                            significant_reactants += 1

                    if significant_reactants >= 2:
                        fragment_couplings += 1
                        print(f"Fragment coupling found at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Consider it a convergent strategy if there are at least 2 fragment couplings
    return fragment_couplings >= 2
