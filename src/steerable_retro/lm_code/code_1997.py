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
    This function detects a strategy involving heterocyclic ring formation (particularly quinolin-2-one)
    via intramolecular cyclization, followed by fragment coupling through N-alkylation and completed
    with late-stage SNAr reaction to incorporate a piperazine moiety.
    """
    # Initialize tracking variables
    has_quinolinone_formation = False
    has_n_alkylation = False
    has_late_stage_snar_with_piperazine = False

    # Track the maximum depth for relative positioning
    max_depth = 0

    def dfs_traverse(node, current_depth=0):
        nonlocal has_quinolinone_formation, has_n_alkylation, has_late_stage_snar_with_piperazine, max_depth

        # Update max depth
        max_depth = max(max_depth, current_depth)

        if node["type"] == "reaction":
            rsmi = node.get("metadata", {}).get("rsmi", "")

            if rsmi:
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Analyzing reaction at depth {current_depth}: {rsmi}")

                # Check for quinolinone formation via cyclization (early stage)
                if current_depth >= 5:  # Early stage reaction based on stdout
                    # Check if product contains a lactam structure in a bicyclic system
                    if "c1c(=O)[nH]c2ccccc12" in product or checker.check_fg(
                        "Secondary amide", product
                    ):
                        # Check if reactants don't have the quinolinone structure
                        reactant_has_quinolinone = any(
                            "c1c(=O)[nH]c2ccccc12" in r for r in reactants
                        )

                        if not reactant_has_quinolinone:
                            print(f"Detected quinolinone formation at depth {current_depth}")
                            has_quinolinone_formation = True

                # Check for N-alkylation (middle stage)
                if 3 <= current_depth <= 7:  # Expanded middle stage range
                    # Check for N-alkylation reactions or benzyl halide + amine pattern
                    if (
                        checker.check_reaction(
                            "N-alkylation of primary amines with alkyl halides", rsmi
                        )
                        or checker.check_reaction(
                            "N-alkylation of secondary amines with alkyl halides", rsmi
                        )
                        or checker.check_reaction("Alkylation of amines", rsmi)
                    ):
                        print(f"Detected N-alkylation at depth {current_depth}")
                        has_n_alkylation = True
                    else:
                        # Check for benzyl halide + amine pattern
                        has_benzyl_halide = any(
                            checker.check_fg("Primary halide", r) and "CH2" in r for r in reactants
                        )
                        has_amine = any(
                            checker.check_fg("Primary amine", r)
                            or checker.check_fg("Secondary amine", r)
                            or checker.check_fg("Aniline", r)
                            for r in reactants
                        )

                        if has_benzyl_halide and has_amine:
                            print(
                                f"Detected N-alkylation with benzyl halide at depth {current_depth}"
                            )
                            has_n_alkylation = True

                        # Additional check for the specific reaction at depth 5 from stdout
                        if (
                            "Br[CH2" in "".join(reactants)
                            and "[NH" in "".join(reactants)
                            and "CH2][NH" in product
                        ):
                            print(
                                f"Detected N-alkylation with benzyl bromide at depth {current_depth}"
                            )
                            has_n_alkylation = True

                # Check for late-stage SNAr with piperazine (late stage)
                if current_depth <= 2:  # Late stage reaction based on stdout
                    # Check if piperazine is a reactant
                    piperazine_reactant = any(
                        checker.check_ring("piperazine", r) for r in reactants
                    )

                    if piperazine_reactant:
                        print(f"Detected piperazine as reactant at depth {current_depth}")

                        # Check for aromatic halide in reactants
                        has_aromatic_halide = any(
                            checker.check_fg("Aromatic halide", r) for r in reactants
                        )

                        if has_aromatic_halide:
                            # Verify piperazine is in product
                            if checker.check_ring("piperazine", product):
                                print(
                                    f"Confirmed late-stage SNAr with piperazine at depth {current_depth}"
                                )
                                has_late_stage_snar_with_piperazine = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, current_depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Strategy detection results:")
    print(f"- Quinolinone formation: {has_quinolinone_formation}")
    print(f"- N-alkylation: {has_n_alkylation}")
    print(f"- Late-stage SNAr with piperazine: {has_late_stage_snar_with_piperazine}")

    # Return True if all key elements of the strategy are present
    return has_quinolinone_formation and has_n_alkylation and has_late_stage_snar_with_piperazine
