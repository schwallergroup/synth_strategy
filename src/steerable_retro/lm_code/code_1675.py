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
    Detects Suzuki coupling reactions in the synthesis route.
    Looks for reactions involving boronic acid and aryl halide.
    """
    result = False

    def dfs_traverse(node, depth=0):
        nonlocal result

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0]

            # Parse reactants
            reactants = [r for r in reactants_smiles.split(".") if r]
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]

            if not all(reactant_mols):
                print("Warning: Could not parse all reactants in reaction")
                return

            # Check for boronic acid pattern
            boronic_acid_pattern = Chem.MolFromSmarts("[#6]-B(O)(O)")
            has_boronic_acid = any(
                mol.HasSubstructMatch(boronic_acid_pattern) for mol in reactant_mols
            )

            # Check for aryl halide pattern
            aryl_halide_pattern = Chem.MolFromSmarts("[#6]:[#6]-[Cl,Br,I]")
            has_aryl_halide = any(
                mol.HasSubstructMatch(aryl_halide_pattern) for mol in reactant_mols
            )

            if has_boronic_acid and has_aryl_halide:
                print(f"Found Suzuki coupling at depth {depth}")
                result = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return result
