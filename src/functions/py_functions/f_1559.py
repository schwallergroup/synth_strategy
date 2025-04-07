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
    Detects a strategy involving fluoropyridine as a key intermediate for SNAr reactions.
    """
    has_fluoropyridine_intermediate = False

    def dfs_traverse(node):
        nonlocal has_fluoropyridine_intermediate

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")

            # Check for fluoropyridine in reactants
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]
            fluoropyridine_pattern = Chem.MolFromSmarts("[c]1[n][c]([F])[c][c][c]1")

            if any(
                mol is not None and mol.HasSubstructMatch(fluoropyridine_pattern)
                for mol in reactant_mols
            ):
                has_fluoropyridine_intermediate = True
                print(f"Found fluoropyridine intermediate in reaction: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(
        f"Fluoropyridine intermediate strategy detected: {has_fluoropyridine_intermediate}"
    )
    return has_fluoropyridine_intermediate
