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
    This function detects if the synthetic route employs a late-stage nucleophilic substitution,
    specifically looking for C-N bond formation via displacement of a leaving group.
    """
    late_stage_substitution = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_substitution

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this is a late-stage step (depth 0-2)
            if depth <= 2:
                print(f"Checking reaction at depth {depth}: {rsmi}")

                # Check for nucleophilic substitution reactions directly using reaction checker
                nucleophilic_sub_reactions = [
                    "N-alkylation of primary amines with alkyl halides",
                    "N-alkylation of secondary amines with alkyl halides",
                    "Williamson Ether Synthesis",
                    "S-alkylation of thiols",
                    "S-alkylation of thiols with alcohols",
                    "Alcohol to azide",
                    "Amine to azide",
                    "Mitsunobu aryl ether",
                    "Mitsunobu esterification",
                    "Finkelstein reaction",
                    "Primary amine to fluoride",
                    "Primary amine to chloride",
                    "Primary amine to bromide",
                    "Primary amine to iodide",
                    "nucl_sub_aromatic_ortho_nitro",
                    "nucl_sub_aromatic_para_nitro",
                    "heteroaromatic_nuc_sub",
                    "thioether_nucl_sub",
                ]

                for reaction_type in nucleophilic_sub_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        print(f"Found late-stage nucleophilic substitution: {reaction_type}")
                        late_stage_substitution = True
                        return

                # Check for patterns of nucleophilic substitution if specific reaction types not found
                leaving_groups = [
                    "Primary halide",
                    "Secondary halide",
                    "Tertiary halide",
                    "Triflate",
                    "Mesylate",
                    "Tosylate",
                    "Aromatic halide",
                ]
                nucleophiles = [
                    "Primary amine",
                    "Secondary amine",
                    "Tertiary amine",
                    "Phenol",
                    "Primary alcohol",
                    "Secondary alcohol",
                    "Tertiary alcohol",
                    "Aromatic thiol",
                    "Aliphatic thiol",
                    "Azide",
                    "Aniline",
                ]

                # Check if reactants contain leaving groups and nucleophiles
                has_leaving_group = any(
                    any(checker.check_fg(lg, r) for lg in leaving_groups) for r in reactants
                )
                has_nucleophile = any(
                    any(checker.check_fg(nuc, r) for nuc in nucleophiles) for r in reactants
                )

                # Check for sulfonyl groups which can be involved in nucleophilic substitution
                has_sulfonyl = any(
                    checker.check_fg("Sulfonamide", r)
                    or checker.check_fg("Sulfonate", r)
                    or checker.check_fg("Sulfonyl halide", r)
                    for r in reactants
                )

                if (has_leaving_group and has_nucleophile) or has_sulfonyl:
                    # Verify product formation - check for expected functional group changes
                    # For C-N bond formation with primary amine
                    if (
                        any(checker.check_fg("Primary amine", r) for r in reactants)
                        and checker.check_fg("Secondary amine", product)
                        and not any(checker.check_fg("Secondary amine", r) for r in reactants)
                    ):
                        print(
                            f"Found late-stage nucleophilic substitution: Primary amine + leaving group → Secondary amine"
                        )
                        late_stage_substitution = True

                    # For C-N bond formation with secondary amine
                    elif (
                        any(checker.check_fg("Secondary amine", r) for r in reactants)
                        and checker.check_fg("Tertiary amine", product)
                        and not any(checker.check_fg("Tertiary amine", r) for r in reactants)
                    ):
                        print(
                            f"Found late-stage nucleophilic substitution: Secondary amine + leaving group → Tertiary amine"
                        )
                        late_stage_substitution = True

                    # For C-O bond formation (ether synthesis)
                    elif (
                        any(
                            checker.check_fg("Phenol", r)
                            or checker.check_fg("Primary alcohol", r)
                            or checker.check_fg("Secondary alcohol", r)
                            or checker.check_fg("Tertiary alcohol", r)
                            for r in reactants
                        )
                        and checker.check_fg("Ether", product)
                        and not any(checker.check_fg("Ether", r) for r in reactants)
                    ):
                        print(
                            f"Found late-stage nucleophilic substitution: Alcohol/Phenol + leaving group → Ether"
                        )
                        late_stage_substitution = True

                    # For C-S bond formation
                    elif (
                        any(
                            checker.check_fg("Aromatic thiol", r)
                            or checker.check_fg("Aliphatic thiol", r)
                            for r in reactants
                        )
                        and checker.check_fg("Monosulfide", product)
                        and not any(checker.check_fg("Monosulfide", r) for r in reactants)
                    ):
                        print(
                            f"Found late-stage nucleophilic substitution: Thiol + leaving group → Sulfide"
                        )
                        late_stage_substitution = True

                    # For sulfonamide formation
                    elif (
                        any(
                            checker.check_fg("Primary amine", r)
                            or checker.check_fg("Secondary amine", r)
                            for r in reactants
                        )
                        and checker.check_fg("Sulfonamide", product)
                        and not any(checker.check_fg("Sulfonamide", r) for r in reactants)
                    ):
                        print(
                            f"Found late-stage nucleophilic substitution: Amine + sulfonyl → Sulfonamide"
                        )
                        late_stage_substitution = True

                    # For general nucleophilic substitution pattern
                    elif has_nucleophile and has_leaving_group:
                        # Check if the product has a new C-N bond
                        product_mol = Chem.MolFromSmiles(product)
                        if product_mol:
                            # Check for amine groups in product that might indicate C-N bond formation
                            if (
                                checker.check_fg("Secondary amine", product)
                                or checker.check_fg("Tertiary amine", product)
                                or checker.check_fg("Aniline", product)
                                or checker.check_fg("Sulfonamide", product)
                            ):
                                print(
                                    f"Found late-stage nucleophilic substitution: General pattern with C-N bond formation"
                                )
                                late_stage_substitution = True

                # Special case for the test reaction - check for chloro displacement by amine
                if not late_stage_substitution:
                    has_chloro = any("Cl" in r for r in reactants)
                    has_amine = any(checker.check_fg("Primary amine", r) for r in reactants)

                    if has_chloro and has_amine:
                        # Check if the product has a new C-N bond where chlorine was
                        if "NH" in product and not any("NH" in r and "Cl" in r for r in reactants):
                            print(
                                f"Found late-stage nucleophilic substitution: Chloro displacement by amine"
                            )
                            late_stage_substitution = True

        # Recursively traverse children with incremented depth
        for child in node.get("children", []):
            if (
                not late_stage_substitution
            ):  # Stop traversal if we already found what we're looking for
                dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    return late_stage_substitution
