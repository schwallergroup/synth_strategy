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
    This function detects a Suzuki coupling reaction involving a nitro-containing aromatic compound.
    """
    has_suzuki_with_nitro = False

    def dfs_traverse(node):
        nonlocal has_suzuki_with_nitro

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Look for boronic acid and halide patterns typical in Suzuki couplings
            boronic_acid_pattern = Chem.MolFromSmarts("B(O)(O)")
            halide_pattern = Chem.MolFromSmarts("[c,C][Br,I,Cl]")
            nitro_pattern = Chem.MolFromSmarts("[N+](=O)[O-]")

            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]

            # Check if reactants contain boronic acid, halide, and nitro group
            has_boronic_acid = any(
                mol and mol.HasSubstructMatch(boronic_acid_pattern) for mol in reactant_mols
            )
            has_halide = any(mol and mol.HasSubstructMatch(halide_pattern) for mol in reactant_mols)
            has_nitro = any(mol and mol.HasSubstructMatch(nitro_pattern) for mol in reactant_mols)

            if has_boronic_acid and has_halide and has_nitro:
                # This is likely a Suzuki coupling with a nitro-containing compound
                has_suzuki_with_nitro = True
                print("Detected Suzuki coupling with nitro-containing compound")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return has_suzuki_with_nitro
