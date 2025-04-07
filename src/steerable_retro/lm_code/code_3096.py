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
    Detects if the synthetic route employs a late-stage Suzuki coupling strategy
    where an aryl halide is coupled with a boronic acid to introduce diversity
    in the penultimate step of the synthesis.
    """
    # Track if we found a Suzuki coupling
    found_suzuki = False
    # Track the depth of the Suzuki coupling
    suzuki_depth = None
    # Track the total depth of the synthesis
    max_depth = 0

    def is_suzuki_coupling(reaction_smiles):
        """Check if a reaction is a Suzuki coupling"""
        # Split into reactants and product
        parts = reaction_smiles.split(">")
        if len(parts) < 3:
            return False

        reactants = parts[0].split(".")
        product = parts[2]

        # Check for aryl halide in reactants
        aryl_halide_pattern = Chem.MolFromSmarts("[c]-[#53,#35,#17]")

        # Check for boronic acid in reactants
        boronic_acid_pattern = Chem.MolFromSmarts("[#5]([#8])[#8]")

        has_aryl_halide = False
        has_boronic_acid = False

        for reactant in reactants:
            try:
                mol = Chem.MolFromSmiles(reactant)
                if mol:
                    if mol.HasSubstructMatch(aryl_halide_pattern):
                        has_aryl_halide = True
                    if mol.HasSubstructMatch(boronic_acid_pattern):
                        has_boronic_acid = True
            except:
                continue

        # If both patterns are found in reactants, it's likely a Suzuki coupling
        return has_aryl_halide and has_boronic_acid

    def dfs_traverse(node, depth=0):
        nonlocal found_suzuki, suzuki_depth, max_depth

        # Update max depth
        max_depth = max(max_depth, depth)

        # Check if this is a reaction node
        if node.get("type") == "reaction":
            # Get reaction SMILES from metadata
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]

                # Check if this is a Suzuki coupling
                if is_suzuki_coupling(rsmi):
                    found_suzuki = True
                    suzuki_depth = depth
                    print(f"Found Suzuki coupling at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Determine if this is a late-stage Suzuki coupling
    # Late stage means it occurs in the first half of the synthesis (lower depth values)
    if found_suzuki and suzuki_depth is not None:
        # If Suzuki occurs in the first half of the synthesis (considering depth)
        is_late_stage = suzuki_depth <= max_depth / 2
        print(
            f"Suzuki depth: {suzuki_depth}, Max depth: {max_depth}, Is late stage: {is_late_stage}"
        )
        return is_late_stage

    return False
