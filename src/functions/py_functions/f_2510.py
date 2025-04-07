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
    This function detects a strategy involving late-stage amide formation.
    It looks for an amide formation reaction in the second half of the synthesis.
    """
    late_amide_formation = False
    max_depth = 0

    # First pass to determine the maximum depth
    def find_max_depth(node, current_depth=0):
        nonlocal max_depth
        max_depth = max(max_depth, current_depth)

        for child in node.get("children", []):
            find_max_depth(child, current_depth + 1)

    # Second pass to check for late-stage amide formation
    def dfs_traverse(node, depth=0):
        nonlocal late_amide_formation

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for amide formation pattern
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
            product_mol = Chem.MolFromSmiles(product) if product else None

            if product_mol:
                # Pattern for acid chloride
                acid_chloride_pattern = Chem.MolFromSmarts("[C](=[O])[Cl]")

                # Pattern for amine
                amine_pattern = Chem.MolFromSmarts("[N;!$(N=*);!$(N#*)]")

                # Pattern for amide in product
                amide_pattern = Chem.MolFromSmarts("[C](=[O])[N;!$(N=*);!$(N#*)]")

                # Check if reactants contain acid chloride and amine
                has_acid_chloride = any(
                    mol and mol.HasSubstructMatch(acid_chloride_pattern)
                    for mol in reactant_mols
                )
                has_amine = any(
                    mol and mol.HasSubstructMatch(amine_pattern)
                    for mol in reactant_mols
                )

                # Check if product contains amide
                has_amide_product = product_mol.HasSubstructMatch(amide_pattern)

                # Check if this is in the second half of the synthesis
                is_late_stage = depth <= (max_depth / 2)

                if (
                    has_acid_chloride
                    and has_amine
                    and has_amide_product
                    and is_late_stage
                ):
                    late_amide_formation = True
                    print(
                        f"Late-stage amide formation detected at depth {depth}: {rsmi}"
                    )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Find maximum depth first
    find_max_depth(route)

    # Then check for late-stage amide formation
    dfs_traverse(route)

    return late_amide_formation
