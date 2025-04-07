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
    This function detects a late-stage Suzuki coupling strategy where an aryl bromide
    and boronic acid form a biaryl system in the second half of the synthesis.
    """
    suzuki_found = False
    suzuki_depth = -1
    max_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal suzuki_found, suzuki_depth, max_depth

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            # Check if this is a Suzuki coupling
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Look for aryl bromide and boronic acid in reactants
                has_aryl_bromide = False
                has_boronic_acid = False

                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol:
                        # Check for aryl bromide
                        if mol.HasSubstructMatch(Chem.MolFromSmarts("[c][Br]")):
                            has_aryl_bromide = True
                        # Check for boronic acid
                        if mol.HasSubstructMatch(Chem.MolFromSmarts("[c][B]([O])[O]")):
                            has_boronic_acid = True

                # Check if product has biaryl
                prod_mol = Chem.MolFromSmiles(product)
                has_biaryl = False
                if prod_mol and prod_mol.HasSubstructMatch(Chem.MolFromSmarts("[c]!@[c]")):
                    has_biaryl = True

                if has_aryl_bromide and has_boronic_acid and has_biaryl:
                    suzuki_found = True
                    suzuki_depth = depth
                    print(f"Found Suzuki coupling at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Consider it late-stage if it occurs in the first half of the synthesis (lower depth)
    is_late_stage = suzuki_found and suzuki_depth <= max_depth / 2

    if is_late_stage:
        print(
            f"Detected late-stage Suzuki coupling strategy at depth {suzuki_depth} (max depth: {max_depth})"
        )

    return is_late_stage
