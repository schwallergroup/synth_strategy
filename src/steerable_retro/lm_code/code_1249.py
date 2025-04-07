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
    Detects if the synthesis involves assembly of multiple heterocyclic structures
    (pyridine, benzodioxole, piperazine).
    """
    heterocycles_found = set()

    def dfs_traverse(node):
        if node["type"] == "mol" and "smiles" in node:
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Check for pyridine
                pyridine_pattern = Chem.MolFromSmarts("c1ccccn1")
                if mol.HasSubstructMatch(pyridine_pattern):
                    heterocycles_found.add("pyridine")

                # Check for benzodioxole
                benzodioxole_pattern = Chem.MolFromSmarts("C1OCOc2ccccc12")
                if mol.HasSubstructMatch(benzodioxole_pattern):
                    heterocycles_found.add("benzodioxole")

                # Check for piperazine
                piperazine_pattern = Chem.MolFromSmarts("N1CCNCC1")
                if mol.HasSubstructMatch(piperazine_pattern):
                    heterocycles_found.add("piperazine")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(f"Found heterocycles: {heterocycles_found}")
    # Return True if at least 2 different heterocycles are found
    return len(heterocycles_found) >= 2
