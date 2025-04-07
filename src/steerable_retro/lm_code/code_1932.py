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
    Detects a 'build-cyclize-couple' strategy where a heterocyclic core is constructed
    in the middle of the synthesis, with cross-coupling reactions used for both initial
    scaffold assembly and late-stage diversification.
    """
    # Track key events
    early_coupling = False
    mid_cyclization = False
    late_coupling = False

    early_coupling_depth = -1
    cyclization_depth = -1
    late_coupling_depth = -1

    # Define SMARTS patterns
    boronic_acid_pattern = Chem.MolFromSmarts("[#5;X3]")
    halide_pattern = Chem.MolFromSmarts("[#53,#35,#9,#17]")

    def is_cross_coupling(reactants, product):
        """Check if a reaction is likely a cross-coupling"""
        # Look for boronic acid in reactants
        has_boronic = any(
            reactant and Chem.MolFromSmiles(reactant).HasSubstructMatch(boronic_acid_pattern)
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
        nonlocal early_coupling, mid_cyclization, late_coupling
        nonlocal early_coupling_depth, cyclization_depth, late_coupling_depth

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for cross-coupling
            if is_cross_coupling(reactants, product):
                if depth <= 1:  # Late stage
                    print(f"Found late-stage cross-coupling at depth {depth}")
                    late_coupling = True
                    late_coupling_depth = depth
                elif depth >= 3:  # Early stage
                    print(f"Found early-stage cross-coupling at depth {depth}")
                    early_coupling = True
                    early_coupling_depth = depth

            # Check for cyclization (ring count increase)
            reactant_rings = sum(count_rings(r) for r in reactants if r)
            product_rings = count_rings(product)

            if product_rings > reactant_rings:
                print(
                    f"Found cyclization at depth {depth}: {reactant_rings} -> {product_rings} rings"
                )
                mid_cyclization = True
                cyclization_depth = depth

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if events occurred in the right order
    correct_order = False
    if early_coupling and mid_cyclization and late_coupling:
        if early_coupling_depth > cyclization_depth > late_coupling_depth:
            correct_order = True

    # Check if the strategy is present
    strategy_present = early_coupling and mid_cyclization and late_coupling and correct_order
    print(f"Build-cyclize-couple strategy detected: {strategy_present}")
    return strategy_present
