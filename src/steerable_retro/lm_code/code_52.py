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

root_data = "/home/andres/Documents/steerable_retro/data"

fg_args = {
    "file_path": f"{root_data}/patterns/functional_groups.json",
    "value_field": "pattern",
    "key_field": "name",
}
reaction_class_args = {
    "file_path": f"{root_data}/patterns/smirks.json",
    "value_field": "smirks",
    "key_field": "name",
}
ring_smiles_args = {
    "file_path": f"{root_data}/patterns/chemical_rings_smiles.json",
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
    This function detects a convergent synthesis strategy involving N-alkylation
    to join a piperazine fragment with a functionalized indole scaffold.
    """
    n_alkylation_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal n_alkylation_detected

        if node["type"] == "reaction":
            # Get reaction SMILES
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Analyzing reaction at depth {depth}: {rsmi}")

                # Check if we have multiple reactants (convergent)
                if len(reactants) >= 2:
                    # Check if this is an N-alkylation reaction (expanded list)
                    is_n_alkylation = (
                        checker.check_reaction(
                            "N-alkylation of primary amines with alkyl halides", rsmi
                        )
                        or checker.check_reaction(
                            "N-alkylation of secondary amines with alkyl halides", rsmi
                        )
                        or checker.check_reaction("Alkylation of amines", rsmi)
                        or checker.check_reaction("N-methylation", rsmi)
                        or checker.check_reaction("Williamson Ether Synthesis", rsmi)
                        or checker.check_reaction("Mitsunobu_tetrazole_1", rsmi)
                    )

                    print(f"Is N-alkylation: {is_n_alkylation}")

                    # Even if not classified as N-alkylation, check for characteristic patterns
                    # Look for reactants with piperazine, indole, and halides/leaving groups
                    piperazine_reactants = []
                    indole_reactants = []
                    halide_reactants = []
                    amine_reactants = []

                    for i, r in enumerate(reactants):
                        if checker.check_ring("piperazine", r):
                            piperazine_reactants.append(i)
                            print(f"Reactant {i} contains piperazine: {r}")

                            # Check if this reactant also has an amine group
                            if (
                                checker.check_fg("Primary amine", r)
                                or checker.check_fg("Secondary amine", r)
                                or checker.check_fg("Tertiary amine", r)
                            ):
                                amine_reactants.append(i)
                                print(f"Reactant {i} contains amine: {r}")

                        if checker.check_ring("indole", r):
                            indole_reactants.append(i)
                            print(f"Reactant {i} contains indole: {r}")

                        if (
                            checker.check_fg("Primary halide", r)
                            or checker.check_fg("Secondary halide", r)
                            or checker.check_fg("Tertiary halide", r)
                            or checker.check_fg("Aromatic halide", r)
                        ):
                            halide_reactants.append(i)
                            print(f"Reactant {i} contains halide: {r}")

                    # Check if product contains both fragments
                    product_has_piperazine = checker.check_ring("piperazine", product)
                    product_has_indole = checker.check_ring("indole", product)

                    print(f"Product has piperazine: {product_has_piperazine}")
                    print(f"Product has indole: {product_has_indole}")

                    # Direct pattern check for N-alkylation
                    # If one reactant has a halide and another has a piperazine or indole with an amine
                    # and the product has both piperazine and indole, it's likely an N-alkylation

                    # Case 1: Standard N-alkylation detection
                    if (
                        is_n_alkylation
                        and product_has_piperazine
                        and product_has_indole
                        and (
                            (
                                len(piperazine_reactants) > 0
                                and len(indole_reactants) > 0
                                and set(piperazine_reactants) != set(indole_reactants)
                            )
                            or (
                                len(halide_reactants) > 0
                                and (
                                    (
                                        len(piperazine_reactants) > 0
                                        and set(piperazine_reactants) != set(halide_reactants)
                                    )
                                    or (
                                        len(indole_reactants) > 0
                                        and set(indole_reactants) != set(halide_reactants)
                                    )
                                )
                            )
                        )
                    ):

                        print(f"Convergent N-alkylation detected at depth {depth}")
                        print(f"Reaction SMILES: {rsmi}")
                        n_alkylation_detected = True

                    # Case 2: Pattern-based detection (even if not classified as N-alkylation)
                    elif (
                        not is_n_alkylation
                        and product_has_piperazine
                        and product_has_indole
                        and len(reactants) >= 2
                        and (
                            (
                                len(piperazine_reactants) > 0
                                and len(indole_reactants) > 0
                                and set(piperazine_reactants) != set(indole_reactants)
                            )
                            or (
                                len(halide_reactants) > 0
                                and (
                                    (
                                        len(piperazine_reactants) > 0
                                        and set(piperazine_reactants) != set(halide_reactants)
                                    )
                                    or (
                                        len(indole_reactants) > 0
                                        and set(indole_reactants) != set(halide_reactants)
                                    )
                                )
                            )
                        )
                    ):

                        # Check if the reaction involves joining fragments
                        # This is a heuristic check - if the product has both fragments but reactants have them separately
                        print(f"Pattern-based convergent N-alkylation detected at depth {depth}")
                        print(f"Reaction SMILES: {rsmi}")
                        n_alkylation_detected = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return n_alkylation_detected
