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
    Detects if the synthesis route involves early-stage fragment coupling via ether linkage.
    Early stage means at high depth in the tree (further from the final product).
    """
    ether_formation_found = False
    ether_formation_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal ether_formation_found, ether_formation_depth

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for ether linkage between aromatic rings in product
            ether_pattern = Chem.MolFromSmarts("c[OX2]Cc")
            product_mol = Chem.MolFromSmiles(product)

            # Check if ether linkage is formed in this reaction
            if product_mol and product_mol.HasSubstructMatch(ether_pattern):
                # Check if reactants don't have this pattern
                reactants_have_ether = False
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and reactant_mol.HasSubstructMatch(ether_pattern):
                        reactants_have_ether = True
                        break

                if not reactants_have_ether:
                    ether_formation_found = True
                    ether_formation_depth = max(ether_formation_depth, depth)
                    print(f"Ether formation detected at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Consider it early-stage if it happens at depth 4 or higher
    return ether_formation_found and ether_formation_depth >= 4
