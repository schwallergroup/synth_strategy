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
    This function detects a strategy involving sequential aromatic functionalization
    followed by a late-stage SNAr coupling with a heterocycle.
    """
    # Track key features
    has_acylation = False
    has_halogenation = False
    has_o_methylation = False
    has_amine_formation = False
    has_late_snar = False
    has_pyrimidine_coupling = False

    # Track the depth of reactions
    min_depth_snar = float("inf")
    max_depth_functionalization = -1

    def dfs_traverse(node, depth=0):
        nonlocal has_acylation, has_halogenation, has_o_methylation, has_amine_formation
        nonlocal has_late_snar, has_pyrimidine_coupling, min_depth_snar, max_depth_functionalization

        if node["type"] == "reaction":
            # Extract reactants and product
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                print(f"Depth {depth}, Reaction: {rsmi}")

                # Check for SNAr reaction
                if (
                    checker.check_reaction("heteroaromatic_nuc_sub", rsmi)
                    or checker.check_reaction("nucl_sub_aromatic_ortho_nitro", rsmi)
                    or checker.check_reaction("nucl_sub_aromatic_para_nitro", rsmi)
                    or (
                        any(checker.check_ring("pyrimidine", r) for r in reactants_smiles)
                        and checker.check_fg("Secondary amine", product_smiles)
                    )
                ):
                    has_late_snar = True
                    min_depth_snar = min(min_depth_snar, depth)
                    print(f"Detected SNAr coupling at depth {depth}")

                    # Check for pyrimidine in reactants
                    for r_smiles in reactants_smiles:
                        if checker.check_ring("pyrimidine", r_smiles):
                            has_pyrimidine_coupling = True
                            print("Detected pyrimidine in SNAr coupling")
                            break

                # Check for amine formation
                if (
                    checker.check_reaction("reductive amination with aldehyde", rsmi)
                    or checker.check_reaction("reductive amination with ketone", rsmi)
                    or checker.check_reaction("reductive amination with alcohol", rsmi)
                    or (
                        checker.check_fg("Secondary amine", product_smiles)
                        and not any(
                            checker.check_fg("Secondary amine", r) for r in reactants_smiles
                        )
                    )
                ):
                    has_amine_formation = True
                    max_depth_functionalization = max(max_depth_functionalization, depth)
                    print(f"Detected amine formation at depth {depth}")

                # Check for O-methylation
                if (
                    checker.check_reaction("O-methylation", rsmi)
                    or checker.check_reaction("Methylation of OH with DMS", rsmi)
                    or (
                        checker.check_fg("Ether", product_smiles)
                        and any(checker.check_fg("Phenol", r) for r in reactants_smiles)
                    )
                ):
                    has_o_methylation = True
                    max_depth_functionalization = max(max_depth_functionalization, depth)
                    print(f"Detected O-methylation at depth {depth}")

                # Check for halogenation
                if (
                    checker.check_reaction("Aromatic bromination", rsmi)
                    or checker.check_reaction("Aromatic chlorination", rsmi)
                    or checker.check_reaction("Aromatic fluorination", rsmi)
                    or checker.check_reaction("Aromatic iodination", rsmi)
                    or (
                        any(
                            checker.check_fg("Aromatic halide", product_smiles)
                            and not checker.check_fg("Aromatic halide", r)
                            for r in reactants_smiles
                        )
                    )
                ):
                    has_halogenation = True
                    max_depth_functionalization = max(max_depth_functionalization, depth)
                    print(f"Detected halogenation at depth {depth}")

                # Check for acylation (Friedel-Crafts)
                if checker.check_reaction("Friedel-Crafts acylation", rsmi) or (
                    checker.check_fg("Ketone", product_smiles)
                    and not any(checker.check_fg("Ketone", r) for r in reactants_smiles)
                ):
                    has_acylation = True
                    max_depth_functionalization = max(max_depth_functionalization, depth)
                    print(f"Detected acylation at depth {depth}")

            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if the strategy is present - SNAr should be at a lower depth (later stage)
    # than at least one functionalization step
    strategy_present = (
        has_late_snar
        and (has_acylation or has_halogenation or has_o_methylation or has_amine_formation)
        and min_depth_snar < max_depth_functionalization
    )

    # Print detailed results for debugging
    print(f"Acylation: {has_acylation}")
    print(f"Halogenation: {has_halogenation}")
    print(f"O-methylation: {has_o_methylation}")
    print(f"Amine formation: {has_amine_formation}")
    print(f"Late SNAr: {has_late_snar}")
    print(f"Pyrimidine coupling: {has_pyrimidine_coupling}")
    print(f"Min depth SNAr: {min_depth_snar}")
    print(f"Max depth functionalization: {max_depth_functionalization}")

    if strategy_present:
        print(
            "Detected strategy: Sequential aromatic functionalization with late-stage SNAr coupling"
        )
    else:
        print("Strategy not fully detected")

    return strategy_present
