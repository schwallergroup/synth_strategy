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
    This function detects convergent synthesis with fragment coupling.
    It looks for reactions where multiple complex fragments are joined.
    """
    fragment_coupling_detected = False

    def dfs_traverse(node):
        nonlocal fragment_coupling_detected

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # Check if we have multiple complex reactants (fragments)
                if len(reactants) >= 2:
                    complex_fragments = 0
                    for reactant in reactants:
                        mol = Chem.MolFromSmiles(reactant)
                        if (
                            mol and mol.GetNumAtoms() > 10
                        ):  # Consider fragments with >10 atoms as complex
                            complex_fragments += 1

                    if complex_fragments >= 2:
                        print(f"Fragment coupling detected: {rsmi}")
                        fragment_coupling_detected = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return fragment_coupling_detected
