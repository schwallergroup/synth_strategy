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
    Detects if the synthetic route employs a convergent synthesis approach
    (multiple fragments joined in a single step)
    """
    convergent_synthesis_detected = False

    def dfs_traverse(node):
        nonlocal convergent_synthesis_detected

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")

            # If there are multiple reactants (more than 1), it might be a convergent step
            if len(reactants) > 1:
                # Check if reactants are substantial fragments (not just reagents)
                substantial_fragments = 0
                for r in reactants:
                    mol = Chem.MolFromSmiles(r) if r else None
                    if mol and mol.GetNumHeavyAtoms() > 5:  # Arbitrary threshold for "substantial"
                        substantial_fragments += 1

                if substantial_fragments >= 2:
                    print(
                        f"Convergent synthesis detected with {substantial_fragments} substantial fragments"
                    )
                    convergent_synthesis_detected = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return convergent_synthesis_detected
