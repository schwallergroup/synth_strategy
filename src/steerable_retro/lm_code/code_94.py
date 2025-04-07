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
    This function detects if cyclopropyl ring formation occurs early in the synthesis.
    """
    early_ring_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal early_ring_formation

        if node["type"] == "reaction" and depth >= 4:  # Early in synthesis (high depth)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0]
                product = rsmi.split(">")[-1]

                # Check for cyclopropyl formation
                if "Cl" in reactants and "Br" in reactants:
                    try:
                        prod_mol = Chem.MolFromSmiles(product)
                        cyclopropyl_patt = Chem.MolFromSmarts("[C]1[C][C]1")
                        if prod_mol and prod_mol.HasSubstructMatch(cyclopropyl_patt):
                            print("Found early cyclopropyl formation at depth", depth)
                            early_ring_formation = True
                    except:
                        pass

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    return early_ring_formation
