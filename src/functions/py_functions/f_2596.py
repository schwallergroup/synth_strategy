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
    This function detects late-stage sulfonamide formation.
    Late stage is defined as occurring in the first half of the synthesis (low depth).
    """
    sulfonamide_formation = False
    max_depth = 0
    sulfonamide_depth = -1

    # First pass to find max depth
    def find_max_depth(node, current_depth=0):
        nonlocal max_depth
        max_depth = max(max_depth, current_depth)

        for child in node.get("children", []):
            find_max_depth(child, current_depth + 1)

    # Second pass to find sulfonamide formation
    def dfs_traverse(node, current_depth=0):
        nonlocal sulfonamide_formation, sulfonamide_depth

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for sulfonyl chloride in reactants
                sulfonyl_chloride = Chem.MolFromSmarts("[S](=O)(=O)Cl")
                amine_pattern = Chem.MolFromSmarts("[NH2]")
                sulfonamide_pattern = Chem.MolFromSmarts("[NH][S](=O)(=O)[#6]")

                has_sulfonyl_chloride = False
                has_amine = False

                for r in reactants:
                    try:
                        mol = Chem.MolFromSmiles(r)
                        if mol:
                            if mol.HasSubstructMatch(sulfonyl_chloride):
                                has_sulfonyl_chloride = True
                            if mol.HasSubstructMatch(amine_pattern):
                                has_amine = True
                    except:
                        continue

                # Check for sulfonamide in product
                try:
                    prod_mol = Chem.MolFromSmiles(product)
                    if prod_mol and prod_mol.HasSubstructMatch(sulfonamide_pattern):
                        if has_sulfonyl_chloride and has_amine:
                            print(
                                f"Detected sulfonamide formation at depth {current_depth}"
                            )
                            sulfonamide_formation = True
                            sulfonamide_depth = current_depth
                except:
                    pass

        for child in node.get("children", []):
            dfs_traverse(child, current_depth + 1)

    find_max_depth(route)
    dfs_traverse(route)

    # Check if sulfonamide formation is in the first half of the synthesis
    if sulfonamide_formation and sulfonamide_depth <= max_depth / 2:
        print(
            f"Sulfonamide formation is late-stage (depth {sulfonamide_depth} out of {max_depth})"
        )
        return True
    return False
