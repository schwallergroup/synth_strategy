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
    This function detects if the synthetic route involves late-stage heterocycle formation.
    Late stage is defined as occurring in the first half of the synthesis (lower depth values).
    """
    max_depth = 0
    heterocycle_formation_depth = None

    # First pass to find max depth
    def find_max_depth(node, current_depth=0):
        nonlocal max_depth
        max_depth = max(max_depth, current_depth)

        for child in node.get("children", []):
            find_max_depth(child, current_depth + 1)

    # Second pass to find heterocycle formations
    def find_heterocycle_formation(node, depth=0):
        nonlocal heterocycle_formation_depth

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Count rings in reactants and product
            product_mol = Chem.MolFromSmiles(product_smiles)
            product_ring_count = 0
            if product_mol:
                product_ring_count = product_mol.GetRingInfo().NumRings()

            reactants_ring_count = 0
            for reactant_smiles in reactants_smiles:
                reactant_mol = Chem.MolFromSmiles(reactant_smiles)
                if reactant_mol:
                    reactants_ring_count += reactant_mol.GetRingInfo().NumRings()

            # Check if a new ring was formed
            if product_ring_count > reactants_ring_count:
                # Check if the new ring contains a heteroatom
                # This is a simplified check - in practice, you'd need to identify which ring is new
                heterocycle_pattern = Chem.MolFromSmarts(
                    "[#6]1[#6,#7,#8,#16][#6,#7,#8,#16][#6,#7,#8,#16][#6,#7,#8,#16][#6]1"
                )
                if product_mol and product_mol.HasSubstructMatch(heterocycle_pattern):
                    heterocycle_formation_depth = depth
                    print(f"Detected heterocycle formation at depth {depth}: {rsmi}")

        for child in node.get("children", []):
            find_heterocycle_formation(child, depth + 1)

    # Run the traversals
    find_max_depth(route)
    find_heterocycle_formation(route)

    # Check if heterocycle formation occurred in the first half of synthesis
    if heterocycle_formation_depth is not None:
        is_late_stage = heterocycle_formation_depth <= (max_depth / 2)
        print(
            f"Heterocycle formation at depth {heterocycle_formation_depth} out of max depth {max_depth}"
        )
        print(f"Is late stage: {is_late_stage}")
        return is_late_stage

    print("No heterocycle formation detected")
    return False
