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
    This function detects if the synthesis follows a linear strategy (each step builds on a single precursor).
    """
    is_linear = True

    def dfs_traverse(node):
        nonlocal is_linear

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # Count complex reactants (more than 10 atoms)
                complex_reactants = 0
                for r in reactants:
                    mol = Chem.MolFromSmiles(r)
                    if mol is not None and mol.GetNumAtoms() > 10:
                        complex_reactants += 1

                # If more than one complex reactant, it's not a linear synthesis
                if complex_reactants > 1:
                    print(
                        f"Detected non-linear step with {complex_reactants} complex reactants"
                    )
                    is_linear = False

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    if is_linear:
        print("Detected linear synthesis strategy")

    return is_linear
