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
    This function detects if an aromatic ether linkage is formed early and maintained throughout the synthesis.
    """
    max_depth = 0
    ether_formation_depth = None

    def dfs_traverse(node, current_depth=0):
        nonlocal max_depth, ether_formation_depth

        # Update max depth
        max_depth = max(max_depth, current_depth)

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Create RDKit mol objects
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
            product_mol = Chem.MolFromSmiles(product) if product else None

            if all(reactant_mols) and product_mol:
                # Check for aromatic ether formation
                aromatic_ether_pattern = Chem.MolFromSmarts("c[#8]c")

                # Check if reactants don't have aromatic ether but product does
                reactants_with_ether = any(
                    mol.HasSubstructMatch(aromatic_ether_pattern) for mol in reactant_mols
                )
                product_has_ether = product_mol.HasSubstructMatch(aromatic_ether_pattern)

                if product_has_ether and not reactants_with_ether:
                    print(f"Detected aromatic ether formation at depth {current_depth}")
                    ether_formation_depth = current_depth

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, current_depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if ether formation occurs early (high depth) and is maintained
    if ether_formation_depth is not None and max_depth > 0:
        # Early formation means high depth value (closer to starting materials)
        if ether_formation_depth >= max_depth * 0.7:  # Consider "early" as in the last 30% of steps
            print(
                f"Confirmed early aromatic ether formation at depth {ether_formation_depth} (max depth: {max_depth})"
            )
            return True

    return False
