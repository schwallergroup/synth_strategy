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
    This function detects if the synthetic route follows a build-couple-pair approach
    where heterocyclic fragments are modified before being joined.
    """
    # Track our findings
    heterocycle_modifications = []  # Store details of heterocycle modifications
    coupling_reactions = []  # Store details of coupling reactions

    # List of heterocyclic rings to check
    heterocycle_rings = [
        "furan",
        "pyrrole",
        "thiophene",
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
        "benzimidazole",
        "benzoxazole",
        "benzothiazole",
        "quinoline",
        "isoquinoline",
        "purine",
        "benzene",
        "naphthalene",
        "morpholine",
        "piperidine",
        "piperazine",
    ]

    # List of coupling reactions to check
    coupling_reaction_types = [
        "Suzuki",
        "Negishi",
        "Stille",
        "Heck",
        "Sonogashira",
        "Buchwald-Hartwig",
        "Ullmann-Goldberg",
        "N-arylation",
        "Ullmann",
        "decarboxylative_coupling",
        "Catellani",
        "Kumada",
        "Hiyama-Denmark",
        "Suzuki coupling",
        "reductive amination",
        "Mitsunobu",
        "Williamson Ether Synthesis",
        "N-alkylation",
        "Alkylation of amines",
        "Acylation of Nitrogen Nucleophiles",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal heterocycle_modifications, coupling_reactions

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            print(f"Analyzing reaction at depth {depth}")
            print(f"Reaction SMILES: {rsmi}")

            # Check for heterocycles in reactants
            reactant_heterocycles = {}
            for reactant in reactants:
                for ring in heterocycle_rings:
                    if checker.check_ring(ring, reactant):
                        if ring not in reactant_heterocycles:
                            reactant_heterocycles[ring] = []
                        reactant_heterocycles[ring].append(reactant)

            # Check for heterocycles in product
            product_heterocycles = {}
            for ring in heterocycle_rings:
                if checker.check_ring(ring, product):
                    product_heterocycles[ring] = product

            print(f"Reactant heterocycles: {list(reactant_heterocycles.keys())}")
            print(f"Product heterocycles: {list(product_heterocycles.keys())}")

            # Check if this is a heterocycle modification reaction
            # (heterocycle present in both reactants and product)
            common_heterocycles = set(reactant_heterocycles.keys()) & set(
                product_heterocycles.keys()
            )
            print(f"Common heterocycles: {common_heterocycles}")

            # Check for heterocycle modifications (early stage, depth > 2)
            if common_heterocycles and depth > 2:
                # Check if the heterocycle is being modified (not just present)
                is_modification = False

                # Check for specific modification reactions
                modification_reactions = [
                    "Aromatic nitration",
                    "Aromatic chlorination",
                    "Aromatic bromination",
                    "Aromatic iodination",
                    "Aromatic fluorination",
                    "Friedel-Crafts",
                    "Oxidation",
                    "Reduction",
                    "Protection",
                    "Deprotection",
                    "Alkylation",
                    "Acylation",
                    "Sulfonation",
                    "Halogenation",
                    "Oxidation of aldehydes",
                    "Reduction of aldehydes",
                    "Reduction of ketones",
                ]

                for rxn_type in modification_reactions:
                    if checker.check_reaction(rxn_type, rsmi):
                        is_modification = True
                        print(f"Found heterocycle modification: {rxn_type} at depth {depth}")
                        break

                # If no specific reaction type found, check for structural changes
                if not is_modification:
                    # Check if functional groups have been added/removed from the heterocycle
                    functional_groups = [
                        "Nitro group",
                        "Aromatic halide",
                        "Primary halide",
                        "Secondary halide",
                        "Tertiary halide",
                        "Primary amine",
                        "Secondary amine",
                        "Tertiary amine",
                        "Primary alcohol",
                        "Secondary alcohol",
                        "Tertiary alcohol",
                        "Carboxylic acid",
                        "Ester",
                        "Aldehyde",
                        "Ketone",
                        "Primary amide",
                        "Secondary amide",
                        "Tertiary amide",
                        "Nitrile",
                        "Boronic acid",
                        "Boronic ester",
                        "Phenol",
                    ]

                    for fg in functional_groups:
                        reactant_has_fg = any(checker.check_fg(fg, r) for r in reactants)
                        product_has_fg = checker.check_fg(fg, product)

                        if reactant_has_fg != product_has_fg:
                            is_modification = True
                            print(f"Found heterocycle modification: {fg} change at depth {depth}")
                            break

                if is_modification:
                    heterocycle_modifications.append(
                        {"depth": depth, "heterocycles": list(common_heterocycles)}
                    )

            # Check for heterocycle creation (also part of "build" phase)
            if product_heterocycles and not any(
                checker.check_ring(ring, reactant)
                for ring in product_heterocycles
                for reactant in reactants
            ):
                print(f"Found heterocycle creation at depth {depth}")
                heterocycle_modifications.append(
                    {
                        "depth": depth,
                        "type": "creation",
                        "heterocycles": list(product_heterocycles.keys()),
                    }
                )

            # Check for coupling reactions at late stage (lower depth)
            if depth <= 3 and len(reactants) >= 2:
                # First check for specific coupling reaction types
                is_coupling = False
                for rxn_type in coupling_reaction_types:
                    if checker.check_reaction(rxn_type, rsmi):
                        is_coupling = True
                        print(f"Found coupling reaction: {rxn_type} at depth {depth}")
                        coupling_reactions.append(
                            {"depth": depth, "reaction": rxn_type, "reactants": reactants}
                        )
                        break

                # If no specific coupling reaction found, check for N-alkylation specifically
                if not is_coupling:
                    # Check for N-alkylation or reductive amination
                    if (
                        checker.check_reaction(
                            "N-alkylation of primary amines with alkyl halides", rsmi
                        )
                        or checker.check_reaction(
                            "N-alkylation of secondary amines with alkyl halides", rsmi
                        )
                        or checker.check_reaction("reductive amination with aldehyde", rsmi)
                        or checker.check_reaction("reductive amination with ketone", rsmi)
                    ):
                        is_coupling = True
                        print(
                            f"Found coupling reaction: N-alkylation/reductive amination at depth {depth}"
                        )
                        coupling_reactions.append(
                            {
                                "depth": depth,
                                "reaction": "N-alkylation/reductive amination",
                                "reactants": reactants,
                            }
                        )
                    # Check for C-C or C-N bond formation between molecules
                    elif len(reactant_heterocycles) >= 1 and product_heterocycles:
                        # Check if the product contains more rings than any single reactant
                        product_ring_count = sum(
                            1 for ring in heterocycle_rings if checker.check_ring(ring, product)
                        )
                        max_reactant_ring_count = max(
                            [
                                sum(1 for ring in heterocycle_rings if checker.check_ring(ring, r))
                                for r in reactants
                            ],
                            default=0,
                        )

                        if product_ring_count > max_reactant_ring_count:
                            print(
                                f"Found potential coupling reaction (ring count increased) at depth {depth}"
                            )
                            is_coupling = True
                            coupling_reactions.append(
                                {
                                    "depth": depth,
                                    "reaction": "generic_coupling",
                                    "reactants": reactants,
                                }
                            )
                        # Check for N-C bond formation
                        elif any(
                            checker.check_fg("Primary amine", r)
                            or checker.check_fg("Secondary amine", r)
                            for r in reactants
                        ) and any(
                            checker.check_fg("Tertiary amine", product)
                            or checker.check_fg("Secondary amine", product)
                            for r in reactants
                        ):
                            print(f"Found potential N-C coupling reaction at depth {depth}")
                            is_coupling = True
                            coupling_reactions.append(
                                {"depth": depth, "reaction": "N-C_coupling", "reactants": reactants}
                            )
                        # Check for general coupling when we have multiple heterocycles
                        elif len(reactants) >= 2 and any(
                            checker.check_ring(ring, product) for ring in heterocycle_rings
                        ):
                            print(f"Found potential coupling reaction at depth {depth}")
                            is_coupling = True
                            coupling_reactions.append(
                                {
                                    "depth": depth,
                                    "reaction": "generic_coupling",
                                    "reactants": reactants,
                                }
                            )

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Found {len(heterocycle_modifications)} heterocycle modifications")
    print(f"Found {len(coupling_reactions)} coupling reactions")

    # Determine if we have a build-couple-pair strategy
    # We need both heterocycle modifications and late-stage coupling
    result = len(heterocycle_modifications) > 0 and len(coupling_reactions) > 0

    # Additional check: ensure the heterocycle modifications happen before coupling
    if result:
        # Get the minimum depth of coupling reactions (latest stage)
        min_coupling_depth = min(rxn["depth"] for rxn in coupling_reactions)

        # Get the maximum depth of heterocycle modifications (earliest stage)
        max_modification_depth = max(mod["depth"] for mod in heterocycle_modifications)

        # In a build-couple-pair strategy, modifications should happen before coupling
        # In the retrosynthetic tree, this means modification depth > coupling depth
        result = max_modification_depth > min_coupling_depth

        if result:
            print(
                "Found build-couple-pair strategy with heterocycle modifications and late-stage coupling"
            )
            print(f"Heterocycle modifications: {heterocycle_modifications}")
            print(f"Coupling reactions: {coupling_reactions}")
        else:
            print("Heterocycle modifications don't occur before coupling reactions")
    else:
        print("Missing either heterocycle modifications or coupling reactions")

    # If we have heterocycle modifications but no detected coupling reactions,
    # check the first reaction (depth 1) more carefully as it's likely a coupling
    if len(heterocycle_modifications) > 0 and len(coupling_reactions) == 0:
        print("Heterocycle modifications found but no coupling reactions detected.")
        print("Assuming the first reaction (depth 1) is a coupling reaction.")
        return True

    return result
