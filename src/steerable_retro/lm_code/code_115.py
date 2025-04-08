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
    Detects if the synthetic route uses a chloroalkyl chain as a leaving group
    for nucleophilic substitution.
    """
    has_chloroalkyl = False

    def dfs_traverse(node):
        nonlocal has_chloroalkyl

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")

            # Check for chloroalkyl pattern in reactants
            chloroalkyl_pattern = Chem.MolFromSmarts("[Cl][CH2][CH2]")

            try:
                for reactant in reactants:
                    if reactant:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and mol.HasSubstructMatch(chloroalkyl_pattern):
                            print("Found chloroalkyl leaving group")
                            has_chloroalkyl = True
                            break
            except:
                pass  # Handle parsing errors gracefully

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return has_chloroalkyl
