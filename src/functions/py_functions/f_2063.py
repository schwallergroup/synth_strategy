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
    This function detects if a heterocycle is formed
    in the late stage of the synthesis.
    """
    # List of heterocycles to check
    heterocycles = [
        "pyrazole",
        "imidazole",
        "oxazole",
        "thiazole",
        "triazole",
        "tetrazole",
        "isoxazole",
        "isothiazole",
        "oxadiazole",
        "thiadiazole",
        "benzoxazole",
        "benzothiazole",
        "benzimidazole",
        "indole",
        "pyrrole",
        "furan",
        "thiophene",
        "pyridine",
        "pyrimidine",
        "pyrazine",
        "pyridazine",
        "quinoline",
        "isoquinoline",
        "piperidine",
        "piperazine",
        "morpholine",
        "thiomorpholine",
    ]

    # Dictionary mapping heterocycles to potential reaction types
    heterocycle_reactions = {
        "pyrazole": ["pyrazole", "{pyrazole}", "Pyrazole formation"],
        "imidazole": ["imidazole", "{imidazole}", "{triaryl-imidazole}"],
        "oxazole": ["oxazole"],
        "thiazole": ["thiazole", "{thiazole}"],
        "triazole": [
            "triazole",
            "{1,2,4-triazole_acetohydrazide}",
            "{1,2,4-triazole_carboxylic-acid/ester}",
            "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
            "Huisgen 1,3 dipolar cycloaddition",
            "{Huisgen_Cu-catalyzed_1,4-subst}",
            "{Huisgen_Ru-catalyzed_1,5_subst}",
        ],
        "tetrazole": [
            "tetrazole",
            "{tetrazole_terminal}",
            "{tetrazole_connect_regioisomere_1}",
            "{tetrazole_connect_regioisomere_2}",
            "Azide-nitrile click cycloaddition to tetrazole",
        ],
        "isoxazole": ["isoxazole"],
        "isothiazole": ["isothiazole"],
        "oxadiazole": [
            "oxadiazole",
            "{oxadiazole}",
            "Huisgen 1,3,4-oxadiazoles from COOH and tetrazole",
        ],
        "thiadiazole": ["thiadiazole"],
        "benzoxazole": [
            "benzoxazole",
            "{benzoxazole}",
            "Benzoxazole formation from aldehyde",
            "Benzoxazole formation from acyl halide",
            "Benzoxazole formation from ester/carboxylic acid",
            "Benzoxazole formation (intramolecular)",
            "{benzoxazole_arom-aldehyde}",
            "{benzoxazole_carboxylic-acid}",
        ],
        "benzothiazole": [
            "benzothiazole",
            "{benzothiazole}",
            "Benzothiazole formation from aldehyde",
            "Benzothiazole formation from acyl halide",
            "Benzothiazole formation from ester/carboxylic acid",
        ],
        "benzimidazole": [
            "benzimidazole",
            "{benzimidazole_derivatives_carboxylic-acid/ester}",
            "{benzimidazole_derivatives_aldehyde}",
            "Benzimidazole formation from aldehyde",
            "Benzimidazole formation from acyl halide",
            "Benzimidazole formation from ester/carboxylic acid",
            "Benzimidazole aldehyde",
        ],
        "indole": [
            "indole",
            "{indole}",
            "{Fischer indole}",
            "Paal-Knorr pyrrole synthesis",
            "{piperidine_indole}",
        ],
        "pyrrole": ["pyrrole", "{Paal-Knorr pyrrole}", "Paal-Knorr pyrrole synthesis"],
        "furan": ["furan", "{benzofuran}"],
        "thiophene": ["thiophene", "{benzothiophene}"],
    }

    # Add missing heterocycles to the reactions dictionary
    for hc in heterocycles:
        if hc not in heterocycle_reactions:
            heterocycle_reactions[hc] = [hc]

    # Track formation depths for each heterocycle
    heterocycle_formation_depths = {hc: -1 for hc in heterocycles}
    max_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal max_depth

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for heterocycle formation
                for heterocycle in heterocycles:
                    # Check if heterocycle exists in reactants and product
                    has_heterocycle_in_reactants = any(
                        checker.check_ring(heterocycle, r) for r in reactants
                    )
                    has_heterocycle_in_product = checker.check_ring(
                        heterocycle, product
                    )

                    # Detect heterocycle formation (not in reactants but in product)
                    if not has_heterocycle_in_reactants and has_heterocycle_in_product:
                        print(
                            f"Potential {heterocycle} formation detected at depth {depth}"
                        )
                        print(f"Reaction SMILES: {rsmi}")

                        # Check for specific heterocycle formation reactions
                        reaction_found = False

                        # Check all possible reaction types for this heterocycle
                        if heterocycle in heterocycle_reactions:
                            for reaction_name in heterocycle_reactions[heterocycle]:
                                if checker.check_reaction(reaction_name, rsmi):
                                    print(
                                        f"Confirmed {heterocycle} formation via {reaction_name} at depth {depth}"
                                    )
                                    heterocycle_formation_depths[heterocycle] = depth
                                    reaction_found = True
                                    break

                        # If no specific reaction was found but we still detected heterocycle formation
                        if not reaction_found:
                            # Check for formation of NOS heterocycles (general reaction)
                            if checker.check_reaction(
                                "Formation of NOS Heterocycles", rsmi
                            ):
                                print(
                                    f"Confirmed {heterocycle} formation via Formation of NOS Heterocycles at depth {depth}"
                                )
                                heterocycle_formation_depths[heterocycle] = depth
                                reaction_found = True

                            # If still no reaction found, accept the heterocycle formation based on structure change
                            if not reaction_found:
                                print(
                                    f"Accepting {heterocycle} formation based on structure change at depth {depth}"
                                )
                                heterocycle_formation_depths[heterocycle] = depth

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Print debug information
    print(f"Max depth: {max_depth}")
    for hc, depth in heterocycle_formation_depths.items():
        if depth != -1:
            print(f"{hc} formation depth: {depth}")

    # Check if any heterocycle formation is in the late stage of the synthesis
    late_stage_threshold = max_depth / 2
    print(f"Late stage threshold (max_depth/2): {late_stage_threshold}")

    for heterocycle, depth in heterocycle_formation_depths.items():
        if depth != -1:
            is_late_stage = depth <= late_stage_threshold
            print(
                f"{heterocycle} formation at depth {depth} is late stage: {is_late_stage}"
            )
            if is_late_stage:
                return True

    return False
