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
    Detects if the synthesis follows a convergent approach where two complex fragments
    are joined together in a coupling reaction.
    """
    convergent_synthesis_detected = False

    def dfs_traverse(node):
        nonlocal convergent_synthesis_detected

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                product = rsmi.split(">")[-1]

                # Split reactants
                reactants = reactants_part.split(".")

                # Check if we have at least 2 complex reactants
                if len(reactants) >= 2:
                    complex_reactants = 0
                    for reactant in reactants:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol:
                            # Define "complex" as having at least 10 heavy atoms
                            if mol.GetNumHeavyAtoms() >= 10:
                                complex_reactants += 1

                    if complex_reactants >= 2:
                        convergent_synthesis_detected = True
                        print(
                            f"Convergent synthesis detected: {complex_reactants} complex reactants joined"
                        )

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return convergent_synthesis_detected
