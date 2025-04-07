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
    Detects if the synthesis involves modification of a heterocycle with morpholine.
    """
    morpholine_modification_detected = False

    # List of heterocycles to check
    heterocycle_rings = [
        "pyrrole",
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
        "benzimidazole",
        "benzoxazole",
        "benzothiazole",
        "furan",
        "thiophene",
        "isoxazole",
        "isothiazole",
    ]

    # List of relevant reaction types for morpholine attachment
    relevant_reactions = [
        "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
        "Buchwald-Hartwig",
        "Ullmann-Goldberg Substitution amine",
        "N-arylation_heterocycles",
        "Nucleophilic substitution",
        "heteroaromatic_nuc_sub",
        "nucl_sub_aromatic_ortho_nitro",
        "nucl_sub_aromatic_para_nitro",
        "N-alkylation of primary amines with alkyl halides",
        "N-alkylation of secondary amines with alkyl halides",
        "Alkylation of amines",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal morpholine_modification_detected

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0]
                product_smiles = rsmi.split(">")[-1]

                # Check for morpholine in reactants and products
                reactants = reactants_smiles.split(".")
                has_morpholine_in_reactants = any(
                    checker.check_ring("morpholine", r) for r in reactants if r
                )
                has_morpholine_in_product = checker.check_ring("morpholine", product_smiles)

                # If morpholine is in either reactants or product, proceed with further checks
                if has_morpholine_in_reactants or has_morpholine_in_product:
                    # Check for heterocycles in product
                    heterocycles_in_product = [
                        ring
                        for ring in heterocycle_rings
                        if checker.check_ring(ring, product_smiles)
                    ]

                    # Check for heterocycles in reactants
                    heterocycles_in_reactants = []
                    for reactant in reactants:
                        if reactant:
                            heterocycles_in_reactants.extend(
                                [
                                    ring
                                    for ring in heterocycle_rings
                                    if checker.check_ring(ring, reactant)
                                ]
                            )

                    # Check if this is a relevant reaction type
                    is_relevant_reaction = any(
                        checker.check_reaction(rxn, rsmi) for rxn in relevant_reactions
                    )

                    # Case 1: Morpholine in reactants, heterocycle in product, relevant reaction type
                    # This covers cases where morpholine is being attached to a heterocycle
                    if (
                        has_morpholine_in_reactants
                        and heterocycles_in_product
                        and is_relevant_reaction
                    ):
                        print(f"Morpholine heterocycle modification detected at depth {depth}")
                        print(f"Reaction SMILES: {rsmi}")
                        print(f"Heterocycles in product: {heterocycles_in_product}")
                        morpholine_modification_detected = True

                    # Case 2: Morpholine in product but not in reactants, and heterocycle present
                    # This covers cases where morpholine is being formed or introduced during the reaction
                    elif (
                        has_morpholine_in_product
                        and not has_morpholine_in_reactants
                        and heterocycles_in_product
                    ):
                        print(f"Morpholine addition to synthesis detected at depth {depth}")
                        print(f"Reaction SMILES: {rsmi}")
                        print(f"Heterocycles in product: {heterocycles_in_product}")
                        morpholine_modification_detected = True

                    # Case 3: Common heterocycles in both reactants and products, and morpholine involved
                    # This covers modifications of heterocycles where morpholine is involved
                    elif (
                        has_morpholine_in_reactants
                        and has_morpholine_in_product
                        and heterocycles_in_reactants
                        and heterocycles_in_product
                        and any(h in heterocycles_in_reactants for h in heterocycles_in_product)
                    ):
                        print(f"Morpholine-heterocycle modification detected at depth {depth}")
                        print(f"Reaction SMILES: {rsmi}")
                        print(f"Heterocycles in product: {heterocycles_in_product}")
                        morpholine_modification_detected = True

                    # Case 4: Check for displacement reactions where morpholine replaces a leaving group on a heterocycle
                    elif (
                        has_morpholine_in_product
                        and heterocycles_in_product
                        and heterocycles_in_reactants
                        and any(h in heterocycles_in_reactants for h in heterocycles_in_product)
                        and (
                            checker.check_fg("Primary halide", reactants_smiles)
                            or checker.check_fg("Secondary halide", reactants_smiles)
                            or checker.check_fg("Tertiary halide", reactants_smiles)
                            or checker.check_fg("Aromatic halide", reactants_smiles)
                        )
                    ):
                        print(f"Morpholine displacement on heterocycle detected at depth {depth}")
                        print(f"Reaction SMILES: {rsmi}")
                        print(f"Heterocycles in product: {heterocycles_in_product}")
                        morpholine_modification_detected = True

            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    return morpholine_modification_detected
