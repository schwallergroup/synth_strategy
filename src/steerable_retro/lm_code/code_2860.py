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
    convergent steps (each reaction has only one non-reagent reactant).
    """
    is_linear = True

    def dfs_traverse(node):
        nonlocal is_linear

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # Count significant reactants (excluding small reagents)
                significant_reactants = 0
                for reactant in reactants:
                    # Skip small molecules/reagents (fewer than 10 atoms)
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.GetNumAtoms() > 10:
                        significant_reactants += 1

                # If more than one significant reactant, it's not a linear synthesis
                if significant_reactants > 1:
                    print(
                        f"Found convergent step with {significant_reactants} significant reactants"
                    )
                    is_linear = False

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return is_linear
