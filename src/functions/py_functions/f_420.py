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
    This function detects if the synthesis involves a late-stage sulfonamide formation
    (in the final or penultimate step).
    """
    sulfonamide_formed = False
    depth_of_formation = float("inf")

    def dfs_traverse(node, depth=0):
        nonlocal sulfonamide_formed, depth_of_formation

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check if this reaction forms a sulfonamide
            product_mol = Chem.MolFromSmiles(product_smiles)
            sulfonamide_pattern = Chem.MolFromSmarts("[N][S](=O)(=O)[C]")

            reactant_has_sulfonamide = False
            for reactant in reactants_smiles:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if reactant_mol and reactant_mol.HasSubstructMatch(sulfonamide_pattern):
                    reactant_has_sulfonamide = True
                    break

            if (
                product_mol
                and product_mol.HasSubstructMatch(sulfonamide_pattern)
                and not reactant_has_sulfonamide
            ):
                sulfonamide_formed = True
                depth_of_formation = min(depth, depth_of_formation)
                print(f"Sulfonamide formation detected at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Consider it late-stage if it happens at depth 0 or 1
    is_late_stage = sulfonamide_formed and depth_of_formation <= 1
    print(
        f"Late-stage sulfonamide formation: {is_late_stage} (formed at depth {depth_of_formation})"
    )
    return is_late_stage
