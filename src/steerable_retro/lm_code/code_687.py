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
    Detects a strategy involving early-stage protections and
    late-stage deprotections in the synthetic route.
    """
    protection_steps = []
    deprotection_steps = []
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal protection_steps, deprotection_steps, max_depth
        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
                product = Chem.MolFromSmiles(product_smiles) if product_smiles else None

                if product and all(r is not None for r in reactants):
                    # Check for protection steps
                    # Acetylation
                    acetyl_chloride_pattern = Chem.MolFromSmarts("C(=O)Cl")
                    amine_pattern = Chem.MolFromSmarts("[NH2]")
                    acetamide_pattern = Chem.MolFromSmarts("NC(=O)C")

                    if (
                        any(r.HasSubstructMatch(acetyl_chloride_pattern) for r in reactants)
                        and any(r.HasSubstructMatch(amine_pattern) for r in reactants)
                        and product.HasSubstructMatch(acetamide_pattern)
                    ):
                        protection_steps.append(depth)
                        print(f"Found acetylation protection at depth {depth}")

                    # Esterification
                    carboxylic_acid_pattern = Chem.MolFromSmarts("[C](=O)[OH]")
                    methyl_ester_pattern = Chem.MolFromSmarts("[C](=O)[O][CH3]")
                    tbutyl_ester_pattern = Chem.MolFromSmarts("[C](=O)[O]C(C)(C)C")

                    if (
                        product.HasSubstructMatch(methyl_ester_pattern)
                        or product.HasSubstructMatch(tbutyl_ester_pattern)
                    ) and any(r.HasSubstructMatch(carboxylic_acid_pattern) for r in reactants):
                        protection_steps.append(depth)
                        print(f"Found esterification protection at depth {depth}")

                    # Check for deprotection steps
                    # Boc deprotection
                    boc_pattern = Chem.MolFromSmarts("[NH][C](=O)[O]C(C)(C)C")

                    if any(
                        r.HasSubstructMatch(boc_pattern) for r in reactants
                    ) and product.HasSubstructMatch(amine_pattern):
                        deprotection_steps.append(depth)
                        print(f"Found Boc deprotection at depth {depth}")

                    # Ester hydrolysis
                    if (
                        any(r.HasSubstructMatch(methyl_ester_pattern) for r in reactants)
                        or any(r.HasSubstructMatch(tbutyl_ester_pattern) for r in reactants)
                    ) and product.HasSubstructMatch(carboxylic_acid_pattern):
                        deprotection_steps.append(depth)
                        print(f"Found ester hydrolysis at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Sort steps by depth
    protection_steps.sort(reverse=True)  # Higher depth = earlier in synthesis
    deprotection_steps.sort()  # Lower depth = later in synthesis

    # Check if we have early protections and late deprotections
    has_early_protection = any(depth > max_depth / 2 for depth in protection_steps)
    has_late_deprotection = any(depth <= max_depth / 2 for depth in deprotection_steps)

    print(f"Protection steps: {protection_steps}")
    print(f"Deprotection steps: {deprotection_steps}")
    print(f"Max depth: {max_depth}")
    print(f"Early protection: {has_early_protection}")
    print(f"Late deprotection: {has_late_deprotection}")

    return has_early_protection and has_late_deprotection
