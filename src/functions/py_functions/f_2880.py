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
    This function detects if the synthetic route involves a convergent synthesis approach
    where two complex fragments are coupled.
    """
    convergent_synthesis_found = False

    def dfs_traverse(node):
        nonlocal convergent_synthesis_found

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for convergent synthesis (multiple complex reactants)
                if len(reactants) >= 2:
                    complex_reactants = 0
                    for reactant in reactants:
                        if reactant:
                            mol = Chem.MolFromSmiles(reactant)
                            if mol:
                                # Consider a molecule complex if it has more than 15 atoms
                                if mol.GetNumAtoms() > 15:
                                    complex_reactants += 1

                    if complex_reactants >= 2:
                        print(
                            "Convergent synthesis detected with multiple complex fragments"
                        )
                        convergent_synthesis_found = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return convergent_synthesis_found
