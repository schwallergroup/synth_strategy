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
    This function detects if the synthesis route involves a late-stage urea formation.
    Late stage is defined as occurring in the first half of the synthesis depth.
    """
    max_depth = 0
    urea_formation_depths = []

    def dfs_traverse(node, current_depth=0):
        nonlocal max_depth, urea_formation_depths

        # Update max depth
        max_depth = max(max_depth, current_depth)

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Create RDKit mol objects
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
            product_mol = Chem.MolFromSmiles(product) if product else None

            if all(reactant_mols) and product_mol:
                # Check for urea formation
                urea_pattern = Chem.MolFromSmarts("[#7][#6](=O)[#7]")

                # Check if reactants don't have urea but product does
                reactants_with_urea = any(
                    mol.HasSubstructMatch(urea_pattern) for mol in reactant_mols
                )
                product_has_urea = product_mol.HasSubstructMatch(urea_pattern)

                if product_has_urea and not reactants_with_urea:
                    print(f"Detected urea formation at depth {current_depth}")
                    urea_formation_depths.append(current_depth)

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, current_depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if urea formation occurs in the first half of the synthesis (late stage)
    if urea_formation_depths and max_depth > 0:
        for depth in urea_formation_depths:
            if depth <= max_depth / 2:  # First half of synthesis (lower depth values)
                print(
                    f"Confirmed late-stage urea formation at depth {depth} (max depth: {max_depth})"
                )
                return True

    return False
