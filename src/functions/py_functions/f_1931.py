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
    Detects a synthetic strategy with cross-coupling reactions at both early and late stages,
    with a heterocycle formation in between.
    """
    # Track if we found cross-couplings at different depths
    early_cross_coupling = False
    late_cross_coupling = False
    heterocycle_formation = False

    # Define SMARTS patterns
    boronic_acid_pattern = Chem.MolFromSmarts("[#5;X3]")
    halide_pattern = Chem.MolFromSmarts("[#53,#35,#9,#17]")

    def is_cross_coupling(reactants, product):
        """Check if a reaction is likely a cross-coupling"""
        # Look for boronic acid in reactants
        has_boronic = any(
            reactant
            and Chem.MolFromSmiles(reactant).HasSubstructMatch(boronic_acid_pattern)
            for reactant in reactants
            if reactant
        )

        # Look for halide in reactants
        has_halide = any(
            reactant and Chem.MolFromSmiles(reactant).HasSubstructMatch(halide_pattern)
            for reactant in reactants
            if reactant
        )

        # If both are present, it's likely a cross-coupling
        return has_boronic and has_halide

    def count_rings(smiles):
        """Count the number of rings in a molecule"""
        if not smiles:
            return 0
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return 0
        return mol.GetRingInfo().NumRings()

    def dfs_traverse(node, depth=0):
        nonlocal early_cross_coupling, late_cross_coupling, heterocycle_formation

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for cross-coupling
            if is_cross_coupling(reactants, product):
                if depth <= 1:  # Late stage (depth 0-1)
                    print(f"Found late-stage cross-coupling at depth {depth}")
                    late_cross_coupling = True
                elif depth >= 3:  # Early stage (depth 3+)
                    print(f"Found early-stage cross-coupling at depth {depth}")
                    early_cross_coupling = True

            # Check for heterocycle formation (ring count increase)
            reactant_rings = sum(count_rings(r) for r in reactants if r)
            product_rings = count_rings(product)

            if product_rings > reactant_rings:
                print(
                    f"Found heterocycle formation at depth {depth}: {reactant_rings} -> {product_rings} rings"
                )
                heterocycle_formation = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if the strategy is present
    strategy_present = (
        early_cross_coupling and late_cross_coupling and heterocycle_formation
    )
    print(f"Bookend cross-coupling strategy detected: {strategy_present}")
    return strategy_present
