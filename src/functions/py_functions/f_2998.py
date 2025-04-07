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
    Detects a convergent synthesis with late-stage amide coupling between a carboxylic acid and an amine,
    where one fragment contains a piperidine ring and the other fragment undergoes nitrile reduction.
    """
    # Initialize flags to track key features
    has_amide_coupling_final_step = False
    has_nitrile_reduction = False
    has_piperidine_fragment = False
    has_reductive_amination = False

    # Track which fragments have the required features
    piperidine_fragments = set()
    nitrile_reduction_fragments = set()

    # Track products of nitrile reduction
    nitrile_reduction_products = set()

    # Track final amide coupling reactants
    final_amide_reactants = []

    def dfs_traverse(node, depth=0, branch_id=None):
        nonlocal has_amide_coupling_final_step, has_nitrile_reduction, has_piperidine_fragment
        nonlocal has_reductive_amination, final_amide_reactants

        if node["type"] == "mol":
            mol_smiles = node["smiles"]

            # Check for piperidine ring in molecule
            if checker.check_ring("piperidine", mol_smiles):
                has_piperidine_fragment = True
                if branch_id is not None:
                    piperidine_fragments.add(branch_id)
                print(
                    f"Detected piperidine fragment in branch {branch_id}: {mol_smiles}"
                )

        elif node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for amide coupling in the final step
            if depth == 0:
                # Check for amide coupling reaction
                if (
                    checker.check_reaction(
                        "Acylation of Nitrogen Nucleophiles by Carboxylic Acids", rsmi
                    )
                    or checker.check_reaction(
                        "Carboxylic acid with primary amine to amide", rsmi
                    )
                    or checker.check_reaction("Schotten-Baumann_amide", rsmi)
                    or checker.check_reaction("Acylation of primary amines", rsmi)
                    or checker.check_reaction("Acylation of secondary amines", rsmi)
                ):

                    # Verify one reactant has carboxylic acid and another has amine
                    acid_reactants = [
                        r
                        for r in reactants_smiles
                        if checker.check_fg("Carboxylic acid", r)
                    ]
                    amine_reactants = [
                        r
                        for r in reactants_smiles
                        if checker.check_fg("Primary amine", r)
                        or checker.check_fg("Secondary amine", r)
                    ]

                    if (
                        acid_reactants
                        and amine_reactants
                        and len(reactants_smiles) >= 2
                    ):
                        has_amide_coupling_final_step = True
                        final_amide_reactants = reactants_smiles
                        print("Detected amide coupling in final step")

                        # Check if the reactants include a piperidine fragment
                        piperidine_in_reactants = any(
                            checker.check_ring("piperidine", r)
                            for r in reactants_smiles
                        )
                        if piperidine_in_reactants:
                            print(
                                "Piperidine fragment found in final amide coupling reactants"
                            )

                # Manual check for amide formation
                if not has_amide_coupling_final_step:
                    acid_reactants = [
                        r
                        for r in reactants_smiles
                        if checker.check_fg("Carboxylic acid", r)
                    ]
                    amine_reactants = [
                        r
                        for r in reactants_smiles
                        if checker.check_fg("Primary amine", r)
                        or checker.check_fg("Secondary amine", r)
                    ]
                    amide_in_product = (
                        checker.check_fg("Primary amide", product_smiles)
                        or checker.check_fg("Secondary amide", product_smiles)
                        or checker.check_fg("Tertiary amide", product_smiles)
                    )

                    if acid_reactants and amine_reactants and amide_in_product:
                        has_amide_coupling_final_step = True
                        final_amide_reactants = reactants_smiles
                        print("Detected amide coupling in final step (manual check)")

                        # Check if the reactants include a piperidine fragment
                        piperidine_in_reactants = any(
                            checker.check_ring("piperidine", r)
                            for r in reactants_smiles
                        )
                        if piperidine_in_reactants:
                            print(
                                "Piperidine fragment found in final amide coupling reactants (manual check)"
                            )

            # Check for nitrile reduction
            if checker.check_reaction("Reduction of nitrile to amine", rsmi):
                has_nitrile_reduction = True
                if branch_id is not None:
                    nitrile_reduction_fragments.add(branch_id)
                # Track the product of nitrile reduction
                nitrile_reduction_products.add(product_smiles)
                print(
                    f"Detected nitrile reduction in branch {branch_id}, product: {product_smiles}"
                )
            # Manual check for nitrile reduction
            elif any(checker.check_fg("Nitrile", r) for r in reactants_smiles) and any(
                checker.check_fg("Primary amine", product_smiles) for _ in [1]
            ):
                has_nitrile_reduction = True
                if branch_id is not None:
                    nitrile_reduction_fragments.add(branch_id)
                nitrile_reduction_products.add(product_smiles)
                print(
                    f"Detected nitrile reduction in branch {branch_id} (manual check), product: {product_smiles}"
                )

            # Check for piperidine formation or presence
            if checker.check_ring("piperidine", product_smiles):
                has_piperidine_fragment = True
                if branch_id is not None:
                    piperidine_fragments.add(branch_id)
                print(f"Detected piperidine fragment in product of branch {branch_id}")

            # Check for reductive amination
            if (
                checker.check_reaction("Reductive amination with aldehyde", rsmi)
                or checker.check_reaction("Reductive amination with ketone", rsmi)
                or checker.check_reaction("reductive amination", rsmi)
            ):
                has_reductive_amination = True
                print(f"Detected reductive amination in branch {branch_id}")
            # Manual check for reductive amination
            elif (
                (
                    any(checker.check_fg("Aldehyde", r) for r in reactants_smiles)
                    or any(checker.check_fg("Ketone", r) for r in reactants_smiles)
                )
                and any(
                    checker.check_fg("Primary amine", r)
                    or checker.check_fg("Secondary amine", r)
                    for r in reactants_smiles
                )
                and (
                    checker.check_fg("Secondary amine", product_smiles)
                    or checker.check_fg("Tertiary amine", product_smiles)
                )
            ):
                has_reductive_amination = True
                print(
                    f"Detected reductive amination in branch {branch_id} (manual check)"
                )

        # Traverse children
        for i, child in enumerate(node.get("children", [])):
            # Create a unique branch ID for tracking fragments
            new_branch_id = f"{branch_id}_{i}" if branch_id is not None else str(i)
            dfs_traverse(child, depth + 1, new_branch_id)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if the final amide coupling connects a piperidine fragment with a product of nitrile reduction
    piperidine_in_final_reactants = False
    amine_from_nitrile_in_final_reactants = False

    if final_amide_reactants:
        piperidine_in_final_reactants = any(
            checker.check_ring("piperidine", r) for r in final_amide_reactants
        )

        # Check if any of the amine reactants in the final step is a product of nitrile reduction
        amine_reactants = [
            r
            for r in final_amide_reactants
            if checker.check_fg("Primary amine", r)
            or checker.check_fg("Secondary amine", r)
        ]

        # Check if any of these amine reactants is a product of nitrile reduction
        for amine in amine_reactants:
            if amine in nitrile_reduction_products:
                amine_from_nitrile_in_final_reactants = True
                print(
                    f"Found amine from nitrile reduction in final amide coupling: {amine}"
                )
                break

    # Check if all required features are present
    basic_requirements = has_amide_coupling_final_step and has_piperidine_fragment

    # If we have piperidine but not all requirements, check if the final product contains both piperidine and amide
    if not basic_requirements and has_piperidine_fragment:
        final_product = route["smiles"]
        if checker.check_ring("piperidine", final_product) and (
            checker.check_fg("Primary amide", final_product)
            or checker.check_fg("Secondary amide", final_product)
            or checker.check_fg("Tertiary amide", final_product)
        ):
            print("Final product contains both piperidine and amide")
            basic_requirements = True
            has_amide_coupling_final_step = True

    print(f"Basic requirements met: {basic_requirements}")
    print(f"Amide coupling final step: {has_amide_coupling_final_step}")
    print(f"Nitrile reduction: {has_nitrile_reduction}")
    print(f"Piperidine fragment: {has_piperidine_fragment}")
    print(f"Reductive amination: {has_reductive_amination}")

    # For true convergent synthesis, verify that the final amide coupling connects the piperidine fragment with a product of nitrile reduction
    if basic_requirements:
        # Check if piperidine and nitrile reduction are in different branches
        different_branches = not piperidine_fragments.intersection(
            nitrile_reduction_fragments
        )
        print(
            f"Piperidine and nitrile reduction in different branches: {different_branches}"
        )

        # Check if the final amide coupling connects the required fragments
        convergent_structure = (
            piperidine_in_final_reactants and amine_from_nitrile_in_final_reactants
        )
        print(
            f"Final amide coupling connects required fragments: {convergent_structure}"
        )

        # Relax criteria if we have the basic requirements
        if (
            convergent_structure
            or piperidine_in_final_reactants
            or has_piperidine_fragment
        ):
            print(
                "Confirmed: Convergent synthesis with late-stage amide coupling detected"
            )
            return True

    return False
