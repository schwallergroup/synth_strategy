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
    This function detects if the synthesis involves late-stage sulfonamide formation.
    """
    found_late_stage_sulfonamide = False

    def dfs_traverse(node, depth=0):
        nonlocal found_late_stage_sulfonamide

        if (
            node["type"] == "reaction" and depth <= 1
        ):  # Check only late-stage reactions (depth 0 or 1)
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this forms a sulfonamide
            try:
                product_mol = Chem.MolFromSmiles(product)

                # Check if product contains sulfonamide group
                if product_mol and product_mol.HasSubstructMatch(
                    Chem.MolFromSmarts("[N;!$(NC=O)]-[S;$(S(=O)(=O))]")
                ):
                    # Check if reactants don't have sulfonamide
                    has_sulfonamide_in_reactants = False
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol and reactant_mol.HasSubstructMatch(
                            Chem.MolFromSmarts("[N;!$(NC=O)]-[S;$(S(=O)(=O))]")
                        ):
                            has_sulfonamide_in_reactants = True
                            break

                    if not has_sulfonamide_in_reactants:
                        found_late_stage_sulfonamide = True
                        print(
                            f"Detected late-stage sulfonamide formation at depth {depth}: {rsmi}"
                        )
            except:
                print(f"Error in processing molecules for sulfonamide check: {rsmi}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Call dfs_traverse on the root node
    dfs_traverse(route)

    print(f"Late-stage sulfonamide formation: {found_late_stage_sulfonamide}")
    return found_late_stage_sulfonamide
