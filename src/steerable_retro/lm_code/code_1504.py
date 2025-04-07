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

from steerable_retro.utils import check, fuzzy_dict
from steerable_retro.utils.check import Check

fg_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/functional_groups.json",
    "value_field": "pattern",
    "key_field": "name",
}
reaction_class_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/smirks.json",
    "value_field": "smirks",
    "key_field": "name",
}
ring_smiles_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/chemical_rings_smiles.json",
    "value_field": "smiles",
    "key_field": "name",
}
functional_groups = fuzzy_dict.FuzzyDict.from_json(**fg_args)
reaction_classes = fuzzy_dict.FuzzyDict.from_json(**reaction_class_args)
ring_smiles = fuzzy_dict.FuzzyDict.from_json(**ring_smiles_args)

checker = check.Check(
    fg_dict=functional_groups, reaction_dict=reaction_classes, ring_dict=ring_smiles
)


def main(route):
    """
    Detects a linear synthesis strategy with early bromination followed by
    functional group transformations and ring formation.
    """
    # Track key transformations and their depths
    bromination_depth = None
    acetylation_depth = None
    o_methylation_depth = None
    ring_formation_depth = None

    print("Starting analysis of synthesis route...")

    def dfs_traverse(node, current_depth=0):
        nonlocal bromination_depth, acetylation_depth, o_methylation_depth, ring_formation_depth

        if node["type"] == "reaction":
            # Get metadata
            metadata = node.get("metadata", {})
            depth = metadata.get("depth", current_depth)
            rsmi = metadata.get("rsmi", "")

            if not rsmi:
                print(f"No reaction SMILES found at depth {depth}")
                return

            print(f"Analyzing reaction at depth {depth}: {rsmi}")

            # Extract reactants and product
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            try:
                # Check for bromination
                if checker.check_reaction("Aromatic bromination", rsmi) or (
                    checker.check_fg("Aromatic halide", product_smiles)
                    and "Br" in product_smiles
                    and not any(
                        checker.check_fg("Aromatic halide", r) and "Br" in r
                        for r in reactants_smiles
                    )
                ):
                    print(f"Bromination detected at depth: {depth}")
                    bromination_depth = depth

                # Check for acetylation (NH2 → NHAc)
                if checker.check_reaction("Acylation of primary amines", rsmi) or (
                    any(checker.check_fg("Primary amine", r) for r in reactants_smiles)
                    and checker.check_fg("Primary amide", product_smiles)
                ):
                    print(f"Acetylation detected at depth: {depth}")
                    acetylation_depth = depth

                # Check for O-methylation (OH → OCH3)
                if checker.check_reaction("O-methylation", rsmi) or (
                    any(checker.check_fg("Phenol", r) for r in reactants_smiles)
                    and not checker.check_fg("Phenol", product_smiles)
                    and checker.check_fg("Ether", product_smiles)
                ):
                    print(f"O-Methylation detected at depth: {depth}")
                    o_methylation_depth = depth

                # Check for benzoxazole formation
                if checker.check_reaction("Benzoxazole formation from aldehyde", rsmi) or (
                    checker.check_ring("benzoxazole", product_smiles)
                    and not any(checker.check_ring("benzoxazole", r) for r in reactants_smiles)
                ):
                    print(f"Benzoxazole formation detected at depth: {depth}")
                    ring_formation_depth = depth

            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, current_depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    print(
        f"Transformation depths - Bromination: {bromination_depth}, Acetylation: {acetylation_depth}, O-Methylation: {o_methylation_depth}, Ring Formation: {ring_formation_depth}"
    )

    # Check if bromination is detected and occurs early (higher depth)
    if bromination_depth is not None:
        # Check if bromination occurs before other transformations
        # Note: Higher depth means earlier in the synthesis (retrosynthetic perspective)
        early_bromination = True

        if acetylation_depth is not None and bromination_depth <= acetylation_depth:
            early_bromination = False

        if o_methylation_depth is not None and bromination_depth <= o_methylation_depth:
            early_bromination = False

        if ring_formation_depth is not None and bromination_depth <= ring_formation_depth:
            early_bromination = False

        # We need at least one other transformation after bromination
        has_subsequent_transformations = (
            acetylation_depth is not None
            or o_methylation_depth is not None
            or ring_formation_depth is not None
        )

        if early_bromination and has_subsequent_transformations:
            print("Linear synthesis with early bromination strategy detected")
            return True
        else:
            print("Bromination detected but not in the expected pattern")
    else:
        print("No bromination detected in the synthesis route")

    return False
