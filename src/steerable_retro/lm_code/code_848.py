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
    This function detects if the synthetic route utilizes aromatic amines (anilines)
    as key nucleophiles.
    """
    aniline_used = False

    def dfs_traverse(node):
        nonlocal aniline_used

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")

            # Check for aniline pattern
            aniline_pattern = Chem.MolFromSmarts("[c][NH2]")

            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if reactant_mol and reactant_mol.HasSubstructMatch(aniline_pattern):
                    aniline_used = True
                    print("Found aromatic amine (aniline) utilization")
                    break

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    return aniline_used
