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
    Detects a synthetic strategy involving protection/deprotection steps
    with linear fragment assembly via amide couplings.
    """
    # Initialize counters and flags
    boc_protection_count = 0
    boc_deprotection_count = 0
    amide_coupling_count = 0
    ester_hydrolysis_count = 0
    late_stage_deprotection = False

    def dfs_traverse(node, depth=0):
        nonlocal boc_protection_count, boc_deprotection_count, amide_coupling_count
        nonlocal ester_hydrolysis_count, late_stage_deprotection

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Convert to RDKit molecules
                reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
                product = Chem.MolFromSmiles(product_smiles) if product_smiles else None

                if product and all(r is not None for r in reactants):
                    # Check for Boc deprotection
                    boc_pattern = Chem.MolFromSmarts("[NH][C](=O)[O]C(C)(C)C")
                    amine_pattern = Chem.MolFromSmarts("[NH2]")

                    if any(
                        r.HasSubstructMatch(boc_pattern) for r in reactants
                    ) and product.HasSubstructMatch(amine_pattern):
                        boc_deprotection_count += 1
                        if depth <= 1:  # Late stage (depth 0 or 1)
                            late_stage_deprotection = True
                            print(f"Found late-stage Boc deprotection at depth {depth}")

                    # Check for amide coupling
                    carboxylic_acid_pattern = Chem.MolFromSmarts("[C](=O)[OH]")
                    amine_pattern = Chem.MolFromSmarts("[NH2]")
                    amide_pattern = Chem.MolFromSmarts("[C](=O)[NH]")

                    if (
                        any(
                            r.HasSubstructMatch(carboxylic_acid_pattern)
                            for r in reactants
                        )
                        and any(r.HasSubstructMatch(amine_pattern) for r in reactants)
                        and product.HasSubstructMatch(amide_pattern)
                    ):
                        amide_coupling_count += 1
                        print(f"Found amide coupling at depth {depth}")

                    # Check for tert-butyl ester hydrolysis
                    tbutyl_ester_pattern = Chem.MolFromSmarts("[C](=O)[O]C(C)(C)C")

                    if any(
                        r.HasSubstructMatch(tbutyl_ester_pattern) for r in reactants
                    ) and product.HasSubstructMatch(carboxylic_acid_pattern):
                        ester_hydrolysis_count += 1
                        print(f"Found tert-butyl ester hydrolysis at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Determine if the strategy is present
    has_protection_deprotection = (
        boc_deprotection_count > 0 and ester_hydrolysis_count > 0
    )
    has_amide_couplings = amide_coupling_count >= 2

    print(
        f"Protection/deprotection count: {boc_deprotection_count + ester_hydrolysis_count}"
    )
    print(f"Amide coupling count: {amide_coupling_count}")
    print(f"Late stage deprotection: {late_stage_deprotection}")

    return (
        has_protection_deprotection and has_amide_couplings and late_stage_deprotection
    )
