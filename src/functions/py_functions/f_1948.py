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
    This function detects the use of a chloro-substituted heterocycle as a key intermediate.
    """
    chloro_heterocycle_pattern = Chem.MolFromSmarts(
        "Cl[c]1[n][c][n][c]2[c][c][c][c][c]12"
    )
    chloro_intermediate_used = False

    def dfs_traverse(node, depth=0):
        nonlocal chloro_intermediate_used

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")

                for reactant in reactants_smiles:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(
                        chloro_heterocycle_pattern
                    ):
                        print(
                            f"Chloro-heterocycle intermediate detected at depth {depth}"
                        )
                        chloro_intermediate_used = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return chloro_intermediate_used
