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
    This function detects if the synthesis follows a linear strategy without
    convergent steps involving multiple complex fragments.
    """
    is_linear = True

    def dfs_traverse(node):
        nonlocal is_linear

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # Count complex reactants (more than 10 atoms)
                complex_reactants = 0
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.GetNumAtoms() > 10:
                        complex_reactants += 1

                # If more than one complex reactant, it's likely a convergent step
                if complex_reactants > 1:
                    print(f"Detected convergent step with {complex_reactants} complex reactants")
                    is_linear = False

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from root
    dfs_traverse(route)

    if is_linear:
        print("Detected linear synthesis strategy")

    return is_linear
