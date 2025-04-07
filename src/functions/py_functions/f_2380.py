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
    Detects late-stage nitrile formation from tertiary amine.
    Looks for a reaction in the first half of the synthesis that introduces a nitrile group.
    """
    nitrile_pattern = Chem.MolFromSmarts("[#6]#[N]")
    tertiary_amine_pattern = Chem.MolFromSmarts("[#6]-[N](-[#6])-[#6]")

    found_nitrile_formation = False
    max_depth = 0
    nitrile_formation_depth = None

    # First pass to find max depth
    def find_max_depth(node, current_depth=0):
        nonlocal max_depth
        max_depth = max(max_depth, current_depth)

        for child in node.get("children", []):
            find_max_depth(child, current_depth + 1)

    # Second pass to find nitrile formation
    def dfs_traverse(node, depth=0):
        nonlocal found_nitrile_formation, nitrile_formation_depth

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if product contains nitrile
                product_mol = Chem.MolFromSmiles(product)
                if product_mol and product_mol.HasSubstructMatch(nitrile_pattern):
                    # Check if any reactant contains tertiary amine
                    for reactant in reactants:
                        reactant_mol = Chem.MolFromSmiles(reactant)
                        if reactant_mol and reactant_mol.HasSubstructMatch(
                            tertiary_amine_pattern
                        ):
                            print(
                                f"Found nitrile formation from tertiary amine at depth {depth}"
                            )
                            found_nitrile_formation = True
                            nitrile_formation_depth = depth
                            break

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    find_max_depth(route)
    dfs_traverse(route)

    # Check if nitrile formation occurs in first half of synthesis (late stage)
    if found_nitrile_formation and nitrile_formation_depth is not None:
        return nitrile_formation_depth <= max_depth / 2

    return False
