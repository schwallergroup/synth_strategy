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
    Detects if a Boc protecting group is maintained throughout most of the synthesis.
    """
    boc_presence_by_depth = {}
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal boc_presence_by_depth, max_depth
        max_depth = max(max_depth, depth)

        if node["type"] == "mol" and "smiles" in node:
            try:
                mol_smiles = node["smiles"]
                # Check if Boc group is present
                has_boc = checker.check_fg("Boc", mol_smiles)
                print(
                    f"Depth {depth}: {'Boc found in' if has_boc else 'No Boc found in'} molecule {mol_smiles}"
                )
                boc_presence_by_depth[depth] = has_boc

                # Additional check for Boc protecting an amine (typical use case)
                if has_boc:
                    # Check if molecule also has an amine that could be protected
                    has_amine = (
                        checker.check_fg("Primary amine", mol_smiles)
                        or checker.check_fg("Secondary amine", mol_smiles)
                        or checker.check_fg("Tertiary amine", mol_smiles)
                    )
                    if has_amine:
                        print(f"Depth {depth}: Molecule has both Boc and amine groups")
            except Exception as e:
                print(f"Error processing molecule at depth {depth}: {e}")

        # Check reaction nodes to verify Boc is maintained through reactions
        elif (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if Boc is present in reactants and products
                reactants_have_boc = any(checker.check_fg("Boc", r) for r in reactants)
                product_has_boc = checker.check_fg("Boc", product)

                # Check if this is a Boc protection or deprotection reaction
                is_boc_protection = checker.check_reaction("Boc amine protection", rsmi)
                is_boc_deprotection = checker.check_reaction(
                    "Boc amine deprotection", rsmi
                )

                # Also check for other Boc-related reactions
                is_boc_protection_explicit = checker.check_reaction(
                    "Boc amine protection explicit", rsmi
                )
                is_boc_protection_anhydride = checker.check_reaction(
                    "Boc amine protection with Boc anhydride", rsmi
                )
                is_boc_protection_ethyl = checker.check_reaction(
                    "Boc amine protection (ethyl Boc)", rsmi
                )
                is_boc_protection_secondary = checker.check_reaction(
                    "Boc amine protection of secondary amine", rsmi
                )
                is_boc_protection_primary = checker.check_reaction(
                    "Boc amine protection of primary amine", rsmi
                )

                # Check for tert-butyl deprotection which might affect Boc
                is_tbutyl_deprotection = checker.check_reaction(
                    "Tert-butyl deprotection of amine", rsmi
                )

                if (
                    is_boc_protection
                    or is_boc_protection_explicit
                    or is_boc_protection_anhydride
                    or is_boc_protection_ethyl
                    or is_boc_protection_secondary
                    or is_boc_protection_primary
                ):
                    print(f"Depth {depth}: Boc protection reaction detected")
                    boc_presence_by_depth[depth] = True
                elif is_boc_deprotection or is_tbutyl_deprotection:
                    print(
                        f"Depth {depth}: Boc deprotection reaction detected - this breaks the maintained protection"
                    )
                    boc_presence_by_depth[depth] = False
                else:
                    # For other reactions, check if Boc is maintained
                    if reactants_have_boc and not product_has_boc:
                        print(f"Depth {depth}: Boc lost in reaction")
                        boc_presence_by_depth[depth] = False
                    elif not reactants_have_boc and product_has_boc:
                        print(f"Depth {depth}: Boc introduced in reaction")
                        boc_presence_by_depth[depth] = True
                    else:
                        print(
                            f"Depth {depth}: Boc {'maintained' if product_has_boc else 'not present'} in reaction"
                        )
                        boc_presence_by_depth[depth] = product_has_boc

                # Check for harsh conditions that might affect Boc
                if "reagents" in node.get("metadata", {}) and isinstance(
                    node["metadata"]["reagents"], str
                ):
                    reagents = node["metadata"]["reagents"].lower()
                    if any(
                        term in reagents
                        for term in [
                            "tfa",
                            "hcl",
                            "h2so4",
                            "strong acid",
                            "high temperature",
                        ]
                    ):
                        print(
                            f"Depth {depth}: Harsh conditions detected that might affect Boc"
                        )
                        boc_presence_by_depth[depth] = False
            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if Boc is present in at least 70% of the synthesis steps
    # Give more weight to late-stage steps (lower depth values)
    boc_steps = 0
    total_steps = 0

    for depth, present in boc_presence_by_depth.items():
        # Late-stage steps have lower depth values
        # Give full weight to final product and immediate precursors, half weight to earlier steps
        weight = 1.0 if depth <= 2 else 0.5

        if present:
            boc_steps += weight
        total_steps += weight

    # Calculate the weighted percentage
    if total_steps > 0:
        boc_percentage = boc_steps / total_steps
        print(
            f"Boc protection maintained with weighted percentage: {boc_percentage:.2f}"
        )

        # Check if Boc is present in at least 70% of steps (weighted)
        if boc_percentage >= 0.7:
            # Additional check: verify Boc is present in the final product (depth 0)
            if 0 in boc_presence_by_depth and boc_presence_by_depth[0]:
                print(
                    "Boc protection maintained throughout synthesis and present in final product"
                )
                return True
            else:
                print("Boc not present in final product, so not maintained throughout")
                return False
        else:
            print(
                f"Boc protection not maintained in enough steps ({boc_percentage:.2f} < 0.7)"
            )
            return False

    print("No steps with Boc protection found")
    return False
