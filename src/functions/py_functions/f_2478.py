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
    This function detects if the final step (or one of the last steps)
    in the synthesis is the formation of a sulfonamide group.
    """
    sulfonamide_at_depth = None
    min_depth = float("inf")

    def dfs_traverse(node, depth=0):
        nonlocal sulfonamide_at_depth, min_depth

        if depth < min_depth:
            min_depth = depth

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for sulfonamide formation
                try:
                    prod_mol = Chem.MolFromSmiles(product)
                    sulfonamide_pattern = Chem.MolFromSmarts("[N]-[S](=[O])(=[O])-[C]")

                    if prod_mol and prod_mol.HasSubstructMatch(sulfonamide_pattern):
                        # Check if reactants don't have sulfonamide
                        has_sulfonamide_in_reactants = False
                        for reactant in reactants:
                            react_mol = Chem.MolFromSmiles(reactant)
                            if react_mol and react_mol.HasSubstructMatch(
                                sulfonamide_pattern
                            ):
                                has_sulfonamide_in_reactants = True
                                break

                        if not has_sulfonamide_in_reactants:
                            print(f"Sulfonamide formation detected at depth {depth}")
                            if (
                                sulfonamide_at_depth is None
                                or depth < sulfonamide_at_depth
                            ):
                                sulfonamide_at_depth = depth
                except:
                    pass

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if sulfonamide formation occurs in the last 30% of steps
    if sulfonamide_at_depth is not None:
        is_late_stage = (
            sulfonamide_at_depth <= min_depth + 1
        )  # Within first two steps (low depth = late stage)
        print(f"Sulfonamide formation is late stage: {is_late_stage}")
        return is_late_stage
    return False
