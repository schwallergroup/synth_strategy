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
    Detects if the route employs a late-stage Suzuki coupling (depth 0 or 1)
    to connect two aromatic fragments.
    """
    suzuki_at_late_stage = False

    def dfs_traverse(node, depth=0):
        nonlocal suzuki_at_late_stage

        if node["type"] == "reaction" and depth <= 1:  # Late stage (depth 0 or 1)
            rsmi = node["metadata"].get("rsmi", "")
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for boronic acid/ester pattern in reactants
            boronic_pattern = Chem.MolFromSmarts("[c]-[B](-[O])-[O]")
            halide_pattern = Chem.MolFromSmarts("[c]-[Cl,Br,I]")

            has_boronic = False
            has_halide = False

            for reactant in reactants:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(boronic_pattern):
                        has_boronic = True
                    if mol and mol.HasSubstructMatch(halide_pattern):
                        has_halide = True
                except:
                    continue

            # Check if product has a new biaryl bond
            if has_boronic and has_halide:
                try:
                    prod_mol = Chem.MolFromSmiles(product)
                    if prod_mol:
                        # This is a simplification - a more robust check would analyze
                        # the actual bond formation between the two aromatic systems
                        biaryl_pattern = Chem.MolFromSmarts("[c]-[c]")
                        if prod_mol.HasSubstructMatch(biaryl_pattern):
                            suzuki_at_late_stage = True
                            print(f"Found late-stage Suzuki coupling at depth {depth}")
                except:
                    pass

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    return suzuki_at_late_stage
