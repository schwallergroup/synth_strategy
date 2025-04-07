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
    This function detects if the synthesis follows a linear strategy (vs convergent).
    """
    is_linear = True
    complex_reactants_count = 0

    def dfs_traverse(node):
        nonlocal is_linear, complex_reactants_count

        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # Count complex reactants (more than 15 atoms)
                complex_count = 0
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.GetNumAtoms() > 15:
                        complex_count += 1

                # If more than one complex reactant, it's likely convergent
                if complex_count > 1:
                    complex_reactants_count += 1
                    print(f"Found step with {complex_count} complex reactants")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # If we found steps with multiple complex reactants, it's not linear
    if complex_reactants_count > 0:
        is_linear = False

    return is_linear
