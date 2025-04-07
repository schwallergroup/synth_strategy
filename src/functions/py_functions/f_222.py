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
    This function detects a strategy involving early-stage diaryl ether formation
    followed by sequential functional group transformations.
    """
    has_diaryl_ether_formation = False
    has_nitrile_hydrolysis = False
    has_amide_formation = False
    has_late_alkylation = False

    def dfs_traverse(node, depth=0):
        nonlocal has_diaryl_ether_formation, has_nitrile_hydrolysis, has_amide_formation, has_late_alkylation

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles]
            product = Chem.MolFromSmiles(product_smiles)

            # Check for diaryl ether formation (early stage, depth >= 3)
            if depth >= 3:
                # Look for formation of diaryl ether (c-O-c)
                diaryl_ether_pattern = Chem.MolFromSmarts("[c]-[O]-[c]")
                if product is not None and product.HasSubstructMatch(
                    diaryl_ether_pattern
                ):
                    # Check if this is a new formation by checking reactants
                    reactant_has_pattern = False
                    for r in reactants:
                        if r is not None and r.HasSubstructMatch(diaryl_ether_pattern):
                            reactant_has_pattern = True
                            break

                    if not reactant_has_pattern:
                        has_diaryl_ether_formation = True
                        print(f"Detected diaryl ether formation at depth {depth}")

            # Check for nitrile hydrolysis
            nitrile_pattern = Chem.MolFromSmarts("[C]#[N]")
            carboxylic_acid_pattern = Chem.MolFromSmarts("[C](=[O])[O]")

            if (
                any(
                    r is not None and r.HasSubstructMatch(nitrile_pattern)
                    for r in reactants
                )
                and product is not None
                and product.HasSubstructMatch(carboxylic_acid_pattern)
            ):
                has_nitrile_hydrolysis = True
                print(f"Detected nitrile hydrolysis at depth {depth}")

            # Check for amide formation
            amide_pattern = Chem.MolFromSmarts("[C](=[O])[N]")
            amine_pattern = Chem.MolFromSmarts("[N;H1,H2]")
            carboxylic_acid_pattern = Chem.MolFromSmarts("[C](=[O])[O]")

            if (
                product is not None
                and product.HasSubstructMatch(amide_pattern)
                and any(
                    r is not None and r.HasSubstructMatch(amine_pattern)
                    for r in reactants
                )
                and any(
                    r is not None and r.HasSubstructMatch(carboxylic_acid_pattern)
                    for r in reactants
                )
            ):
                has_amide_formation = True
                print(f"Detected amide formation at depth {depth}")

            # Check for late-stage alkylations (N-methylation or O-methylation)
            if depth <= 1:  # Late stage
                # N-methylation: secondary amine to tertiary amine
                if product is not None:
                    n_methyl_pattern = Chem.MolFromSmarts("[N]-[CH3]")
                    if product.HasSubstructMatch(n_methyl_pattern) and not any(
                        r is not None and r.HasSubstructMatch(n_methyl_pattern)
                        for r in reactants
                    ):
                        has_late_alkylation = True
                        print(f"Detected N-methylation at depth {depth}")

                # O-methylation: phenol to methoxy
                o_methyl_pattern = Chem.MolFromSmarts("[c]-[O]-[CH3]")
                phenol_pattern = Chem.MolFromSmarts("[c]-[OH]")

                if (
                    product is not None
                    and product.HasSubstructMatch(o_methyl_pattern)
                    and any(
                        r is not None and r.HasSubstructMatch(phenol_pattern)
                        for r in reactants
                    )
                ):
                    has_late_alkylation = True
                    print(f"Detected O-methylation at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Strategy is present if we have diaryl ether formation early and at least two of the other features
    features_count = sum(
        [has_nitrile_hydrolysis, has_amide_formation, has_late_alkylation]
    )
    strategy_present = has_diaryl_ether_formation and features_count >= 2

    print(f"Early diaryl ether formation strategy detected: {strategy_present}")
    return strategy_present
