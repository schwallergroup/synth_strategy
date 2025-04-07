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
    Detects the use of silyl protection strategy for amines.
    """
    silyl_protection_detected = False

    def dfs_traverse(node):
        nonlocal silyl_protection_detected

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")

            # Check for silyl-protected amine pattern
            silyl_amine_pattern = Chem.MolFromSmarts("[#7]-[Si]")

            for reactant in reactants:
                mol = Chem.MolFromSmiles(reactant)
                if mol and mol.HasSubstructMatch(silyl_amine_pattern):
                    print("Silyl protection strategy detected for amine")
                    silyl_protection_detected = True
                    break

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return silyl_protection_detected
