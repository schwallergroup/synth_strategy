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
    This function detects a sequential transformation from nitro to amine to amide.
    """
    # Track transformations by depth
    transformations = {}

    def dfs_traverse(node):
        if node["type"] == "reaction":
            depth = node.get("depth", 0)

            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles]
            product = Chem.MolFromSmiles(product_smiles)

            if product is not None:
                nitro_pattern = Chem.MolFromSmarts("[N+](=O)[O-]")
                amine_pattern = Chem.MolFromSmarts("[NH2]")
                amide_pattern = Chem.MolFromSmarts("[NH][C](=O)")

                # Check for nitro reduction
                if any(
                    r is not None and r.HasSubstructMatch(nitro_pattern)
                    for r in reactants
                ) and product.HasSubstructMatch(amine_pattern):
                    transformations[depth] = "nitro_to_amine"

                # Check for amide formation
                if any(
                    r is not None and r.HasSubstructMatch(amine_pattern)
                    for r in reactants
                ) and product.HasSubstructMatch(amide_pattern):
                    transformations[depth] = "amine_to_amide"

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Check if we have both transformations in the correct sequence
    has_nitro_to_amine = "nitro_to_amine" in transformations.values()
    has_amine_to_amide = "amine_to_amide" in transformations.values()

    # Get depths for each transformation
    nitro_to_amine_depth = next(
        (d for d, t in transformations.items() if t == "nitro_to_amine"), None
    )
    amine_to_amide_depth = next(
        (d for d, t in transformations.items() if t == "amine_to_amide"), None
    )

    # Check if the sequence is correct (nitro→amine→amide)
    correct_sequence = (
        has_nitro_to_amine
        and has_amine_to_amide
        and (
            nitro_to_amine_depth is not None
            and amine_to_amide_depth is not None
            and nitro_to_amine_depth
            > amine_to_amide_depth  # Remember: higher depth = earlier in synthesis
        )
    )

    if correct_sequence:
        print(
            f"Detected nitro→amine→amide sequence: nitro reduction at depth {nitro_to_amine_depth}, amide formation at depth {amine_to_amide_depth}"
        )

    return correct_sequence
