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
    This function detects if the synthesis follows a linear strategy (vs. convergent).
    """
    is_linear = True

    def dfs_traverse(node):
        nonlocal is_linear

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")

                # If more than 2 complex reactants, consider it convergent
                complex_reactants = 0
                for reactant_smiles in reactants_smiles:
                    # Consider a reactant complex if it has more than 15 atoms
                    reactant_mol = Chem.MolFromSmiles(reactant_smiles)
                    if reactant_mol and reactant_mol.GetNumAtoms() > 15:
                        complex_reactants += 1

                if complex_reactants > 2:
                    print(
                        f"Convergent step detected with {complex_reactants} complex reactants: {rsmi}"
                    )
                    is_linear = False
            except Exception as e:
                print(f"Error in linear synthesis detection: {e}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return is_linear
