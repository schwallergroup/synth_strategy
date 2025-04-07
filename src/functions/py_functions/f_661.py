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
    This function detects a Suzuki coupling strategy where at least one of the coupling partners
    contains a heterocyclic system (like indole, pyridine, etc.).
    """
    found_heterocycle_suzuki = False

    def dfs_traverse(node, depth=0):
        nonlocal found_heterocycle_suzuki

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Convert to RDKit molecules
            reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
            product = Chem.MolFromSmiles(product_smiles) if product_smiles else None

            if product and all(reactants):
                # Check for Suzuki coupling
                boronic_pattern = Chem.MolFromSmarts("[B][O]")
                halide_pattern = Chem.MolFromSmarts("[c][Br,I,Cl]")

                has_boronic = any(
                    mol.HasSubstructMatch(boronic_pattern) for mol in reactants
                )
                has_halide = any(
                    mol.HasSubstructMatch(halide_pattern) for mol in reactants
                )

                # Check for heterocycles
                indole_pattern = Chem.MolFromSmarts("[c]1[c][c][c][c]2[c]1[n][c][c]2")
                pyridine_pattern = Chem.MolFromSmarts("[c]1[c][c][c][n][c]1")
                pyrrole_pattern = Chem.MolFromSmarts("[c]1[c][c][n][c]1")

                heterocycle_patterns = [
                    indole_pattern,
                    pyridine_pattern,
                    pyrrole_pattern,
                ]

                has_heterocycle = any(
                    any(
                        mol.HasSubstructMatch(pattern)
                        for pattern in heterocycle_patterns
                    )
                    for mol in reactants
                )

                if (
                    has_boronic and has_halide and has_heterocycle and depth <= 1
                ):  # Late-stage
                    found_heterocycle_suzuki = True
                    print(
                        f"Found heterocycle-containing Suzuki coupling at depth {depth}"
                    )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return found_heterocycle_suzuki
