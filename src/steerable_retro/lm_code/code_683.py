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
    This function detects a synthetic strategy involving nitro reduction followed by
    protection of the resulting amine (e.g., with a carbamate).
    """
    nitro_reduction_detected = False
    amine_protection_detected = False
    nitro_reduction_depth = -1
    amine_protection_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal nitro_reduction_detected, amine_protection_detected, nitro_reduction_depth, amine_protection_depth

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0]
                product_smiles = rsmi.split(">")[-1]

                # Parse reactants and product
                reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles.split(".")]
                product = Chem.MolFromSmiles(product_smiles)

                if None not in reactants and product is not None:
                    # Check for nitro reduction
                    nitro_pattern = Chem.MolFromSmarts("[N+](=O)[O-]")
                    amine_pattern = Chem.MolFromSmarts("[NH2]")

                    if product.HasSubstructMatch(amine_pattern) and any(
                        r.HasSubstructMatch(nitro_pattern) for r in reactants
                    ):
                        nitro_reduction_detected = True
                        nitro_reduction_depth = depth
                        print(f"Nitro reduction detected at depth {depth}")

                    # Check for amine protection (carbamate formation)
                    carbamate_pattern = Chem.MolFromSmarts("[#8]=[#6]([#7])[#8]")
                    acid_chloride_pattern = Chem.MolFromSmarts("[#6](=[#8])[Cl]")

                    if product.HasSubstructMatch(carbamate_pattern) and any(
                        r.HasSubstructMatch(acid_chloride_pattern) for r in reactants
                    ):
                        amine_protection_detected = True
                        amine_protection_depth = depth
                        print(f"Amine protection (carbamate formation) detected at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if nitro reduction occurs before amine protection in the retrosynthetic direction
    # (which means amine protection occurs after nitro reduction in the forward direction)
    return (
        nitro_reduction_detected
        and amine_protection_detected
        and nitro_reduction_depth > amine_protection_depth
    )
