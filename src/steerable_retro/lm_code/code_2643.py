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
    Detects a synthetic strategy involving elaboration of a pyridine core
    while preserving a nitro group throughout the synthesis.
    """
    # Track presence of nitro group through reactions
    reactions_with_nitro = 0
    total_reactions = 0

    # SMARTS patterns
    pyridine_pattern = Chem.MolFromSmarts("[n;r6]1ccccc1")
    nitro_pattern = Chem.MolFromSmarts("[N+](=[O])[O-]")

    def dfs_traverse(node):
        nonlocal reactions_with_nitro, total_reactions

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                total_reactions += 1
                rsmi = node["metadata"]["rsmi"]
                reactants_str = rsmi.split(">")[0]
                product_str = rsmi.split(">")[-1]

                reactants = [Chem.MolFromSmiles(r) for r in reactants_str.split(".") if r]
                product = Chem.MolFromSmiles(product_str) if product_str else None

                if product and any(reactant for reactant in reactants):
                    # Check if both reactant and product have nitro group
                    reactant_has_nitro = any(
                        reactant and reactant.HasSubstructMatch(nitro_pattern)
                        for reactant in reactants
                    )
                    product_has_nitro = product and product.HasSubstructMatch(nitro_pattern)

                    if reactant_has_nitro and product_has_nitro:
                        reactions_with_nitro += 1
                        print(
                            f"Detected reaction preserving nitro group (total: {reactions_with_nitro})"
                        )

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Determine if the strategy is present (nitro preserved in at least 3 reactions)
    strategy_present = reactions_with_nitro >= 3 and reactions_with_nitro == total_reactions

    print(f"Pyridine elaboration with preserved nitro strategy detected: {strategy_present}")
    print(f"- Reactions with nitro preservation: {reactions_with_nitro}/{total_reactions}")

    return strategy_present
