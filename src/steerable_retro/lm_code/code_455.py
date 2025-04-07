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
    Detects if the synthesis uses a convergent approach with coupling of two heterocyclic systems.
    """
    found_convergent_coupling = False

    def dfs_traverse(node):
        nonlocal found_convergent_coupling

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                if len(reactants) >= 2:
                    # Patterns for common heterocycles
                    heterocycle_patterns = [
                        Chem.MolFromSmarts("[#16]1[#6][#7][#6][#6]1"),  # thiazole
                        Chem.MolFromSmarts("[#7]1[#6][#6][#6][#6][#6]1"),  # pyridine
                        Chem.MolFromSmarts("[#7]1[#6][#7][#6][#6]1"),  # imidazole
                        Chem.MolFromSmarts("[#7]1[#6][#6][#7][#6][#6]1"),  # pyrimidine
                        Chem.MolFromSmarts("[#7]1[#6][#6][#6]2[#6][#6][#6][#6][#6]12"),  # indole
                    ]

                    heterocycle_count = 0
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol:
                            for pattern in heterocycle_patterns:
                                if reactant_mol.HasSubstructMatch(pattern):
                                    heterocycle_count += 1
                                    break

                    if heterocycle_count >= 2:
                        print("Found convergent coupling of heterocyclic systems")
                        found_convergent_coupling = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return found_convergent_coupling
