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
    This function detects a strategy involving early amide formation from ester and amine,
    followed by late-stage SNAr reaction to incorporate a nitrogen heterocycle.
    """
    # Initialize tracking variables
    has_early_amide_formation = False
    has_late_stage_snar = False

    # List of nitrogen heterocycles to check
    n_heterocycles = [
        "pyridine",
        "pyrazole",
        "imidazole",
        "oxazole",
        "thiazole",
        "pyrimidine",
        "pyrazine",
        "pyridazine",
        "triazole",
        "tetrazole",
        "indole",
        "quinoline",
        "isoquinoline",
        "purine",
        "benzimidazole",
        "benzoxazole",
        "benzothiazole",
        "indazole",
        "benzotriazole",
        "thiazole",
    ]

    # Store all reactions in order of traversal
    all_reactions = []

    # Calculate depth of a node in the synthesis tree
    def calculate_depth(node, root):
        if node == root:
            return 0

        def find_depth_recursive(current, target, current_depth=0):
            if current == target:
                return current_depth

            for child in current.get("children", []):
                result = find_depth_recursive(child, target, current_depth + 1)
                if result is not None:
                    return result
            return None

        return find_depth_recursive(root, node)

    def dfs_traverse(node):
        nonlocal has_early_amide_formation, has_late_stage_snar

        if node["type"] == "reaction":
            # Try to get depth from metadata or calculate it
            depth = node.get("metadata", {}).get("depth", None)
            if depth is None:
                # For this implementation, we'll use position in all_reactions as a proxy for depth
                depth = len(all_reactions)
            else:
                depth = int(depth)

            rsmi = node.get("metadata", {}).get("rsmi", "")

            if rsmi:
                all_reactions.append((depth, rsmi))
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Analyzing reaction at depth {depth}: {rsmi}")

                # Check for amide formation
                amide_formation_reactions = [
                    "Ester with primary amine to amide",
                    "Ester with secondary amine to amide",
                    "Aminolysis of esters",
                    "Carboxylic acid with primary amine to amide",
                    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                    "Acyl chloride with secondary amine to amide",
                    "Schotten-Baumann_amide",
                    "Acylation of primary amines",
                    "Acylation of secondary amines",
                    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                ]

                is_amide_formation = any(
                    checker.check_reaction(rxn, rsmi) for rxn in amide_formation_reactions
                )

                # Check if product has amide
                has_amide_product = (
                    checker.check_fg("Secondary amide", product)
                    or checker.check_fg("Primary amide", product)
                    or checker.check_fg("Tertiary amide", product)
                )

                # Check if reactants include ester/acid and amine
                has_ester_or_acid = any(
                    checker.check_fg(fg, r)
                    for r in reactants
                    for fg in ["Ester", "Carboxylic acid", "Acyl halide"]
                )
                has_amine = any(
                    checker.check_fg(fg, r)
                    for r in reactants
                    for fg in ["Primary amine", "Secondary amine", "Aniline"]
                )

                # Check if reactants don't have amide but product does (amide formation)
                reactants_have_amide = any(
                    checker.check_fg(fg, r)
                    for r in reactants
                    for fg in ["Primary amide", "Secondary amide", "Tertiary amide"]
                )

                if is_amide_formation or (
                    has_amide_product
                    and has_ester_or_acid
                    and has_amine
                    and not reactants_have_amide
                ):
                    print(f"✓ Detected amide formation at depth {depth}")
                    # Mark as early if depth >= 4 or if it's in the first half of reactions
                    if depth >= 4 or depth >= len(all_reactions) // 2:
                        print(f"✓ This is considered early-stage amide formation")
                        has_early_amide_formation = True

                # Check for SNAr or N-arylation reactions
                snar_reactions = [
                    "heteroaromatic_nuc_sub",
                    "nucl_sub_aromatic_ortho_nitro",
                    "nucl_sub_aromatic_para_nitro",
                    "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
                    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
                    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
                    "Buchwald-Hartwig",
                    "N-arylation_heterocycles",
                    "Goldberg coupling",
                    "Ullmann-Goldberg Substitution amine",
                ]

                is_snar_or_arylation = any(
                    checker.check_reaction(rxn, rsmi) for rxn in snar_reactions
                )

                # Check for nitrogen heterocycles in product
                product_n_heterocycles = [
                    ring for ring in n_heterocycles if checker.check_ring(ring, product)
                ]

                # Check for leaving groups in reactants
                has_leaving_group = any(
                    checker.check_fg(fg, r)
                    for r in reactants
                    for fg in ["Aromatic halide", "Triflate", "Tosylate", "Mesylate"]
                )

                # Check for nitrogen nucleophiles in reactants
                has_nitrogen_nucleophile = any(
                    checker.check_fg(fg, r)
                    for r in reactants
                    for fg in ["Primary amine", "Secondary amine", "Tertiary amine", "Aniline"]
                )

                # Manual check for SNAr-like transformation
                if (
                    not is_snar_or_arylation
                    and has_leaving_group
                    and has_nitrogen_nucleophile
                    and product_n_heterocycles
                ):
                    # Check if a new C-N bond is formed in a heterocycle
                    print(
                        f"Detected potential SNAr-like transformation with N-heterocycle at depth {depth}"
                    )
                    is_snar_or_arylation = True

                if is_snar_or_arylation and product_n_heterocycles:
                    print(
                        f"✓ Detected SNAr or N-arylation at depth {depth} with nitrogen heterocycle: {product_n_heterocycles}"
                    )
                    # Mark as late-stage if depth <= 2 or if it's in the second half of reactions
                    if depth <= 2 or depth < len(all_reactions) // 2:
                        print(f"✓ This is considered late-stage SNAr")
                        has_late_stage_snar = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # If we have reactions but couldn't determine early/late stages properly,
    # try to infer the strategy from the sequence of reactions
    if all_reactions and not (has_early_amide_formation and has_late_stage_snar):
        # Sort reactions by depth/position
        all_reactions.sort(key=lambda x: x[0])

        # Check if amide formation occurs before SNAr/N-arylation
        amide_indices = []
        snar_indices = []

        for i, (depth, rsmi) in enumerate(all_reactions):
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for amide formation
            amide_formation_reactions = [
                "Ester with primary amine to amide",
                "Ester with secondary amine to amide",
                "Aminolysis of esters",
                "Carboxylic acid with primary amine to amide",
                "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                "Acyl chloride with secondary amine to amide",
                "Schotten-Baumann_amide",
                "Acylation of primary amines",
                "Acylation of secondary amines",
                "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
            ]

            is_amide_formation = any(
                checker.check_reaction(rxn, rsmi) for rxn in amide_formation_reactions
            )

            # Check if product has amide
            has_amide_product = (
                checker.check_fg("Secondary amide", product)
                or checker.check_fg("Primary amide", product)
                or checker.check_fg("Tertiary amide", product)
            )

            # Check if reactants include ester/acid and amine
            has_ester_or_acid = any(
                checker.check_fg(fg, r)
                for r in reactants
                for fg in ["Ester", "Carboxylic acid", "Acyl halide"]
            )
            has_amine = any(
                checker.check_fg(fg, r)
                for r in reactants
                for fg in ["Primary amine", "Secondary amine", "Aniline"]
            )

            # Check if reactants don't have amide but product does (amide formation)
            reactants_have_amide = any(
                checker.check_fg(fg, r)
                for r in reactants
                for fg in ["Primary amide", "Secondary amide", "Tertiary amide"]
            )

            if is_amide_formation or (
                has_amide_product and has_ester_or_acid and has_amine and not reactants_have_amide
            ):
                amide_indices.append(i)

            # Check for SNAr or N-arylation
            snar_reactions = [
                "heteroaromatic_nuc_sub",
                "nucl_sub_aromatic_ortho_nitro",
                "nucl_sub_aromatic_para_nitro",
                "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
                "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
                "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
                "Buchwald-Hartwig",
                "N-arylation_heterocycles",
                "Goldberg coupling",
                "Ullmann-Goldberg Substitution amine",
            ]

            is_snar_or_arylation = any(checker.check_reaction(rxn, rsmi) for rxn in snar_reactions)

            # Check for nitrogen heterocycles in product
            product_n_heterocycles = [
                ring for ring in n_heterocycles if checker.check_ring(ring, product)
            ]

            # Check for leaving groups in reactants
            has_leaving_group = any(
                checker.check_fg(fg, r)
                for r in reactants
                for fg in ["Aromatic halide", "Triflate", "Tosylate", "Mesylate"]
            )

            # Check for nitrogen nucleophiles in reactants
            has_nitrogen_nucleophile = any(
                checker.check_fg(fg, r)
                for r in reactants
                for fg in ["Primary amine", "Secondary amine", "Tertiary amine", "Aniline"]
            )

            # Manual check for SNAr-like transformation
            if (
                not is_snar_or_arylation
                and has_leaving_group
                and has_nitrogen_nucleophile
                and product_n_heterocycles
            ):
                is_snar_or_arylation = True

            if is_snar_or_arylation and product_n_heterocycles:
                snar_indices.append(i)

        # Check if we have both reaction types and amide formation comes before SNAr
        if amide_indices and snar_indices:
            min_amide_idx = min(amide_indices)
            max_snar_idx = max(snar_indices)

            if min_amide_idx < max_snar_idx:
                print(
                    f"Strategy detected from reaction sequence: amide formation at position {min_amide_idx} followed by SNAr at position {max_snar_idx}"
                )
                has_early_amide_formation = True
                has_late_stage_snar = True

    # Return True if both key elements of the strategy are present
    result = has_early_amide_formation and has_late_stage_snar
    print(
        f"Strategy detected: early amide formation: {has_early_amide_formation}, late-stage SNAr: {has_late_stage_snar}"
    )
    return result
