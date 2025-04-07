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
    where two complex fragments are joined together.
    """
    convergent_synthesis_found = False

    def dfs_traverse(node, depth=0):
        nonlocal convergent_synthesis_found

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                reactants = reactants_part.split(".")

                # Check if we have multiple complex reactants (convergent synthesis)
                if len(reactants) >= 2:
                    complex_reactants = 0
                    for reactant in reactants:
                        mol = Chem.MolFromSmiles(reactant)
                        if (
                            mol and mol.GetNumAtoms() > 10
                        ):  # Consider reactants with >10 atoms as complex
                            complex_reactants += 1

                    if complex_reactants >= 2:
                        print(
                            f"Found convergent synthesis at depth {depth} with {complex_reactants} complex fragments"
                        )
                        convergent_synthesis_found = True

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return convergent_synthesis_found
