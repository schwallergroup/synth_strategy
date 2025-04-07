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
    This function detects a strategy involving indole ring formation via
    nitro group reduction.
    """
    has_nitro_group = False
    has_indole_formation = False
    nitro_depth = None
    indole_depth = None

    def dfs_traverse(node, depth=0):
        nonlocal has_nitro_group, has_indole_formation, nitro_depth, indole_depth

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for nitro group in reactants
            nitro_pattern = Chem.MolFromSmarts("[#6]-[N+](=[O])[O-]")
            reactants_mols = [Chem.MolFromSmiles(r) for r in reactants]

            if any(
                r and r.HasSubstructMatch(nitro_pattern) for r in reactants_mols if r
            ):
                has_nitro_group = True
                nitro_depth = depth
                print(f"Nitro group detected at depth {depth}")

            # Check for indole formation
            product_mol = Chem.MolFromSmiles(product)
            indole_pattern = Chem.MolFromSmarts("[c]1[c][c][c][c]2[c]1[nH][c][c]2")

            if product_mol and product_mol.HasSubstructMatch(indole_pattern):
                # Check if indole is formed in this reaction (not present in reactants)
                indole_in_reactants = any(
                    r and r.HasSubstructMatch(indole_pattern)
                    for r in reactants_mols
                    if r
                )
                if not indole_in_reactants:
                    has_indole_formation = True
                    indole_depth = depth
                    print(f"Indole formation detected at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if nitro group is present and then indole is formed
    if has_nitro_group and has_indole_formation:
        if nitro_depth is not None and indole_depth is not None:
            if (
                nitro_depth >= indole_depth
            ):  # Same or higher depth means earlier or same step
                print("Strategy detected: Indole formation via nitro reduction")
                return True

    return False
