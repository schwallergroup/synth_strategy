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
    This function detects if the synthesis involves a late-stage ether formation
    (in the final or penultimate step).
    """
    ether_formation_depth = None
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal ether_formation_depth, max_depth

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for ether formation
                product_mol = Chem.MolFromSmiles(product)
                if product_mol:
                    ether_patt = Chem.MolFromSmarts("[#6]-[#8]-[#6]")
                    if product_mol.HasSubstructMatch(ether_patt):
                        # Check if reactants don't have the ether pattern
                        ether_in_reactants = False
                        for reactant in reactants:
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            if reactant_mol and reactant_mol.HasSubstructMatch(
                                ether_patt
                            ):
                                ether_in_reactants = True

                        if not ether_in_reactants:
                            ether_formation_depth = depth
                            print(f"Detected ether formation at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if ether formation is in the final or penultimate step (depth 0 or 1)
    if ether_formation_depth is not None and ether_formation_depth <= 1:
        print(f"Confirmed late-stage ether formation at depth {ether_formation_depth}")
        return True
    return False
