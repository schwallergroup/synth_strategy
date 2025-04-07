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
    Detects if the synthesis incorporates a 4-cyanoaniline fragment.
    """
    cyanoaniline_used = False

    def dfs_traverse(node):
        nonlocal cyanoaniline_used

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # Check for 4-cyanoaniline fragment
                cyanoaniline_pattern = Chem.MolFromSmarts("c1cc(N)ccc1C#N")

                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(cyanoaniline_pattern):
                        cyanoaniline_used = True
                        print("4-Cyanoaniline fragment detected")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    return cyanoaniline_used
