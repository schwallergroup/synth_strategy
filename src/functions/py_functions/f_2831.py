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
    This function detects if the synthesis follows a linear strategy without convergent steps.
    """
    is_linear = True

    def dfs_traverse(node):
        nonlocal is_linear

        if node["type"] == "reaction":
            # Count number of non-trivial reactants (complex molecules)
            complex_reactants = 0

            if node.get("metadata", {}).get("rsmi"):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                for reactant in reactants:
                    # Consider a reactant complex if it has more than 10 atoms
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.GetNumAtoms() > 10:
                        complex_reactants += 1

            # If a reaction has more than one complex reactant, it's likely convergent
            if complex_reactants > 1:
                is_linear = False
                print(
                    f"Detected convergent step in reaction {node.get('metadata', {}).get('ID', '')}"
                )

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return is_linear
