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
    Detects a synthetic strategy involving sequential transformation of
    nitro group to amine to sulfonamide on an aromatic ring.
    """
    # Initialize tracking variables
    nitro_reduction_found = False
    sulfonamide_formation_found = False
    correct_sequence = False
    nitro_reduction_depth = -1
    sulfonamide_formation_depth = -1

    # SMARTS patterns
    nitro_pattern = Chem.MolFromSmarts("[#6]-[N+](=[O])[O-]")
    amine_pattern = Chem.MolFromSmarts("c-[NH2]")
    sulfonamide_pattern = Chem.MolFromSmarts("[#7]-[S](=O)(=O)-[#6]")

    def dfs_traverse(node, depth=0):
        nonlocal nitro_reduction_found, sulfonamide_formation_found
        nonlocal nitro_reduction_depth, sulfonamide_formation_depth

        if node["type"] == "reaction":
            # Extract reactants and products
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Convert to RDKit molecules
            reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
            product = Chem.MolFromSmiles(product_smiles) if product_smiles else None

            if product and reactants:
                # Check for nitro reduction
                if not nitro_reduction_found:
                    reactant_has_nitro = any(
                        r and r.HasSubstructMatch(nitro_pattern) for r in reactants
                    )
                    product_has_amine = product and product.HasSubstructMatch(
                        amine_pattern
                    )
                    product_has_no_nitro = product and not product.HasSubstructMatch(
                        nitro_pattern
                    )

                    if (
                        reactant_has_nitro
                        and product_has_amine
                        and product_has_no_nitro
                    ):
                        nitro_reduction_found = True
                        nitro_reduction_depth = depth
                        print(f"Nitro reduction detected at depth {depth}")

                # Check for sulfonamide formation
                if not sulfonamide_formation_found:
                    reactant_has_amine = any(
                        r and r.HasSubstructMatch(amine_pattern) for r in reactants
                    )
                    product_has_sulfonamide = product and product.HasSubstructMatch(
                        sulfonamide_pattern
                    )

                    if reactant_has_amine and product_has_sulfonamide:
                        sulfonamide_formation_found = True
                        sulfonamide_formation_depth = depth
                        print(f"Sulfonamide formation detected at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if the transformations occurred in the correct sequence
    if nitro_reduction_found and sulfonamide_formation_found:
        # In retrosynthetic direction, nitro reduction should be at higher depth than sulfonamide formation
        correct_sequence = nitro_reduction_depth > sulfonamide_formation_depth

    strategy_present = (
        nitro_reduction_found and sulfonamide_formation_found and correct_sequence
    )
    print(f"Nitro to sulfonamide sequence strategy detected: {strategy_present}")
    return strategy_present
