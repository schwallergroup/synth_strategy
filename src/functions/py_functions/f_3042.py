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
    This function detects a synthetic strategy involving transformation of aromatic halides
    to nitrogen-containing functional groups (amines, nitro) or interconversion between
    different nitrogen-containing functional groups on aromatic rings.
    """
    has_aromatic_n_functionalization = False

    def dfs_traverse(node, depth=0):
        nonlocal has_aromatic_n_functionalization

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            print(f"Analyzing reaction at depth {depth}: {rsmi}")

            # Check for aromatic rings in both reactants and products
            has_aromatic_reactant = any(
                checker.check_ring("benzene", r) for r in reactants
            )
            has_aromatic_product = checker.check_ring("benzene", product)

            if has_aromatic_reactant and has_aromatic_product:
                # Check for nitrogen-containing functional groups
                n_containing_reactants = any(
                    any(
                        checker.check_fg(fg, r)
                        for fg in [
                            "Aniline",
                            "Nitro group",
                            "Azide",
                            "Nitrile",
                            "Primary amide",
                            "Secondary amide",
                            "Tertiary amide",
                            "Urea",
                            "Thiourea",
                        ]
                    )
                    for r in reactants
                )

                n_containing_product = any(
                    checker.check_fg(fg, product)
                    for fg in [
                        "Aniline",
                        "Nitro group",
                        "Azide",
                        "Nitrile",
                        "Primary amide",
                        "Secondary amide",
                        "Tertiary amide",
                        "Urea",
                        "Thiourea",
                    ]
                )

                # Check for aromatic halide in reactants
                has_aromatic_halide_reactant = any(
                    checker.check_fg("Aromatic halide", r) for r in reactants
                )

                # Case 1: Aromatic halide to nitrogen-containing group (forward)
                if has_aromatic_halide_reactant and n_containing_product:
                    if (
                        checker.check_reaction(
                            "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
                            rsmi,
                        )
                        or checker.check_reaction(
                            "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
                            rsmi,
                        )
                        or checker.check_reaction(
                            "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)", rsmi
                        )
                        or checker.check_reaction("N-arylation_heterocycles", rsmi)
                        or checker.check_reaction("Goldberg coupling", rsmi)
                        or checker.check_reaction(
                            "Goldberg coupling aryl amine-aryl chloride", rsmi
                        )
                        or checker.check_reaction(
                            "Goldberg coupling aryl amide-aryl chloride", rsmi
                        )
                        or checker.check_reaction(
                            "Ullmann-Goldberg Substitution amine", rsmi
                        )
                        or checker.check_reaction("Buchwald-Hartwig", rsmi)
                        or checker.check_reaction(
                            "Formation of Azides from halogens", rsmi
                        )
                    ):
                        has_aromatic_n_functionalization = True
                        print(
                            f"Detected aromatic halide to N-containing group conversion (forward): {rsmi}"
                        )

                # Case 2: Nitro reduction to amine (forward)
                nitro_to_amine = any(
                    checker.check_fg("Nitro group", r) for r in reactants
                ) and checker.check_fg("Aniline", product)
                if nitro_to_amine:
                    if checker.check_reaction(
                        "Reduction of nitro groups to amines", rsmi
                    ):
                        has_aromatic_n_functionalization = True
                        print(
                            f"Detected nitro group reduction to amine (forward): {rsmi}"
                        )

                # Case 3: Amine to nitro (retrosynthetic)
                amine_to_nitro = checker.check_fg("Aniline", product) and any(
                    checker.check_fg("Nitro group", r) for r in reactants
                )
                if amine_to_nitro:
                    if checker.check_reaction(
                        "Reduction of nitro groups to amines", rsmi
                    ):
                        has_aromatic_n_functionalization = True
                        print(f"Detected amine to nitro group (retrosynthetic): {rsmi}")

                # Case 4: Aromatic nitration
                if (
                    checker.check_reaction("Aromatic nitration with HNO3", rsmi)
                    or checker.check_reaction("Aromatic nitration with NO3 salt", rsmi)
                    or checker.check_reaction("Aromatic nitration with NO2 salt", rsmi)
                    or checker.check_reaction("Aromatic nitration with alkyl NO2", rsmi)
                ):
                    if checker.check_fg("Nitro group", product):
                        has_aromatic_n_functionalization = True
                        print(f"Detected aromatic nitration: {rsmi}")

                # Case 5: Other nitrogen functionalization reactions
                if n_containing_reactants or n_containing_product:
                    if (
                        checker.check_reaction(
                            "Formation of Azides from boronic acids", rsmi
                        )
                        or checker.check_reaction("Amine to azide", rsmi)
                        or checker.check_reaction(
                            "Reductive amination with aldehyde", rsmi
                        )
                        or checker.check_reaction(
                            "Reductive amination with ketone", rsmi
                        )
                        or checker.check_reaction(
                            "Reductive amination with alcohol", rsmi
                        )
                    ):
                        has_aromatic_n_functionalization = True
                        print(
                            f"Detected nitrogen functionalization on aromatic system: {rsmi}"
                        )

                # Case 6: Check for oxidation/reduction between different N-containing groups
                if n_containing_reactants and n_containing_product:
                    # This covers interconversion between different nitrogen functional groups
                    has_aromatic_n_functionalization = True
                    print(
                        f"Detected interconversion between N-containing groups: {rsmi}"
                    )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return has_aromatic_n_functionalization
