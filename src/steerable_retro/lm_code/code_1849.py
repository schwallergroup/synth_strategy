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
    This function detects a linear functionalization strategy of an indole scaffold
    with multiple sequential modifications.
    """
    indole_scaffold_present = False
    modification_count = 0

    def dfs_traverse(node):
        nonlocal indole_scaffold_present, modification_count

        # Check for indole scaffold in molecules
        if node["type"] == "mol" and "smiles" in node:
            indole_pattern = Chem.MolFromSmarts("c1ccc2[nH]ccc2c1")
            chloro_indole_pattern = Chem.MolFromSmarts("c1cc(Cl)cc2[nH]ccc12")

            try:
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol:
                    if mol.HasSubstructMatch(indole_pattern):
                        indole_scaffold_present = True
            except:
                print("Error in SMILES processing for indole scaffold detection")

        # Count modifications in reactions
        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                modification_count += 1

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    # Linear functionalization strategy requires indole scaffold and multiple modifications
    return indole_scaffold_present and modification_count >= 4
