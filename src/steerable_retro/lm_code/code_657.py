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
    This function detects if the final step (or one of the last steps) is a biaryl coupling.
    """
    biaryl_coupling_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal biaryl_coupling_depth

        if node["type"] == "reaction":
            # Extract reactants and products
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0]
            product = rsmi.split(">")[-1]

            # Check for biaryl coupling patterns
            # Look for boronic acid/ester + aryl halide → biaryl
            boronic_pattern = Chem.MolFromSmarts("[c]-[B]")
            aryl_halide_pattern = Chem.MolFromSmarts("[c]-[Cl,Br,I]")
            biaryl_pattern = Chem.MolFromSmarts("[c]-[c]")

            reactants_mol = Chem.MolFromSmiles(reactants)
            product_mol = Chem.MolFromSmiles(product)

            if reactants_mol and product_mol:
                if (
                    reactants_mol.HasSubstructMatch(boronic_pattern)
                    and reactants_mol.HasSubstructMatch(aryl_halide_pattern)
                    and product_mol.HasSubstructMatch(biaryl_pattern)
                ):
                    biaryl_coupling_depth = depth
                    print(f"Biaryl coupling found at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if biaryl coupling was found at a low depth (late stage)
    if biaryl_coupling_depth >= 0 and biaryl_coupling_depth <= 1:
        print(f"Late-stage biaryl coupling strategy detected at depth {biaryl_coupling_depth}")
        return True
    return False
