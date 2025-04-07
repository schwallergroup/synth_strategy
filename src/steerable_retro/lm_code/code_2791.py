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
    This function detects if the synthetic route involves multiple heterocyclic compounds,
    specifically focusing on pyridine rings.
    """
    heterocycle_count = 0
    unique_heterocycles = set()

    def dfs_traverse(node, depth=0):
        nonlocal heterocycle_count, unique_heterocycles

        if node["type"] == "mol" and "smiles" in node:
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Pyridine pattern
                pyridine_pattern = Chem.MolFromSmarts("c1ccncc1")

                # Check for pyridine substructure
                if mol.HasSubstructMatch(pyridine_pattern):
                    # Use a simplified representation of the molecule to avoid counting the same core multiple times
                    # This is a simplification - in a real implementation you might want a more sophisticated approach
                    mol_hash = Chem.MolToSmiles(mol, isomericSmiles=False)[:20]
                    if mol_hash not in unique_heterocycles:
                        unique_heterocycles.add(mol_hash)
                        heterocycle_count += 1
                        print(f"Found heterocyclic compound at depth {depth}")

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Return True if at least 2 unique heterocyclic compounds were detected
    return heterocycle_count >= 2
