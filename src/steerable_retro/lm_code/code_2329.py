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
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")

                # Count complex reactants (non-reagents)
                complex_reactants = 0
                for reactant in reactants_smiles:
                    # Simple heuristic: complex molecules have more than 15 atoms
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.GetNumAtoms() > 15:
                        complex_reactants += 1

                # If more than one complex reactant, it's likely a convergent step
                if complex_reactants > 1:
                    is_linear = False
                    print("Convergent synthesis step detected")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from root
    dfs_traverse(route)

    return is_linear
