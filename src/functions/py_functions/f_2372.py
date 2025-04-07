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
    by checking if multiple fragments are combined in a single reaction.
    """
    convergent_synthesis = False

    def dfs_traverse(node):
        nonlocal convergent_synthesis

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")

            # Check if we have multiple significant reactants (convergent)
            if len(reactants) > 1:
                significant_reactants = 0
                for reactant in reactants:
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if (
                            mol and mol.GetNumHeavyAtoms() > 6
                        ):  # Consider only substantial fragments
                            significant_reactants += 1
                    except:
                        continue

                if significant_reactants >= 2:
                    print(
                        "Convergent synthesis detected with multiple significant reactants in RSMI:",
                        rsmi,
                    )
                    convergent_synthesis = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return convergent_synthesis
