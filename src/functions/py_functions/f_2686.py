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
    This function detects early-stage spirocyclic system formation involving
    an indole and a small ring.
    """
    spirocycle_formed = False
    indole_pattern = Chem.MolFromSmarts(
        "[#6]1:[#6]:[#6]:[#6]2:[#7]:[#6]:[#6]:[#6]:[#6]:2:[#6]:1"
    )

    def dfs_traverse(node, depth=0):
        nonlocal spirocycle_formed

        if node["type"] == "reaction" and depth >= 3:  # Early stage (high depth)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if product contains a spirocyclic system
                product_mol = Chem.MolFromSmiles(product)
                if product_mol:
                    # Check for spirocyclic system
                    for atom in product_mol.GetAtoms():
                        if (
                            atom.IsInRingSize(3)
                            and len([r for r in atom.GetNeighbors() if r.IsInRing()])
                            >= 3
                        ):
                            # Check if indole is present
                            if product_mol.HasSubstructMatch(indole_pattern):
                                print(f"Detected spirocycle formation at depth {depth}")
                                spirocycle_formed = True
                                break

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return spirocycle_formed
