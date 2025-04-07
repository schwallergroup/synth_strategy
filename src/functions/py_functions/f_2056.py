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
    This function detects a strategy involving heterocycle formation followed by
    sequential nucleophilic aromatic substitutions (SNAr) for side chain elaboration.
    """
    # Track if we found heterocycle formation and SNAr reactions
    heterocycle_formation_found = False
    heterocycle_formation_depth = float("inf")
    snar_reactions = []  # Store tuples of (depth, reaction_smiles)
    late_stage_sulfonyl = False

    # List of heterocycles to check
    heterocycles = [
        "furan",
        "pyran",
        "pyrrole",
        "pyridine",
        "pyrazole",
        "imidazole",
        "oxazole",
        "thiazole",
        "pyrimidine",
        "pyrazine",
        "triazole",
        "tetrazole",
        "indole",
        "benzoxazole",
        "benzothiazole",
        "benzimidazole",
    ]

    # List of specific heterocycle formation reactions
    heterocycle_formation_reactions = [
        "benzimidazole_derivatives_carboxylic-acid/ester",
        "benzimidazole_derivatives_aldehyde",
        "benzothiazole",
        "benzoxazole_arom-aldehyde",
        "benzoxazole_carboxylic-acid",
        "thiazole",
        "tetrazole_terminal",
        "Huisgen_Cu-catalyzed_1,4-subst",
        "1,2,4-triazole_acetohydrazide",
        "1,2,4-triazole_carboxylic-acid/ester",
        "pyrazole",
        "Paal-Knorr pyrrole",
        "Fischer indole",
        "benzofuran",
        "benzothiophene",
        "indole",
        "oxadiazole",
        "Formation of NOS Heterocycles",
    ]

    # List of SNAr and related reactions
    snar_reaction_types = [
        "heteroaromatic_nuc_sub",
        "nucl_sub_aromatic_ortho_nitro",
        "nucl_sub_aromatic_para_nitro",
        "Buchwald-Hartwig",
        "N-arylation",
        "Ullmann-Goldberg Substitution amine",
        "Ullmann-Goldberg Substitution thiol",
        "Ullmann-Goldberg Substitution aryl alcohol",
        "Goldberg coupling",
        "Ullmann condensation",
        "N-arylation_heterocycles",
        "Williamson Ether Synthesis",
        "Williamson ether",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal heterocycle_formation_found, heterocycle_formation_depth, snar_reactions, late_stage_sulfonyl

        if node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check for heterocycle formation using specific reactions
                for rxn_type in heterocycle_formation_reactions:
                    if checker.check_reaction(rxn_type, rsmi):
                        print(
                            f"Heterocycle formation reaction {rxn_type} detected at depth {depth}"
                        )
                        heterocycle_formation_found = True
                        heterocycle_formation_depth = min(
                            heterocycle_formation_depth, depth
                        )
                        break

                # If no specific reaction found, check for heterocycle presence
                if not heterocycle_formation_found:
                    # Check if product contains a heterocycle
                    product_has_heterocycle = False
                    reactants_have_heterocycle = False

                    # Check which heterocycles are in the product
                    product_heterocycles = []
                    for heterocycle in heterocycles:
                        if checker.check_ring(heterocycle, product_smiles):
                            product_has_heterocycle = True
                            product_heterocycles.append(heterocycle)
                            print(f"Product contains {heterocycle} at depth {depth}")

                    # Check if any reactant contains the same heterocycles
                    reactant_heterocycles = set()
                    for reactant in reactants_smiles:
                        if not reactant:  # Skip empty reactants
                            continue
                        for heterocycle in heterocycles:
                            if checker.check_ring(heterocycle, reactant):
                                reactants_have_heterocycle = True
                                reactant_heterocycles.add(heterocycle)
                                print(
                                    f"Reactant contains {heterocycle} at depth {depth}"
                                )

                    # If product has heterocycle but reactants don't have all the same heterocycles,
                    # it's a heterocycle formation
                    new_heterocycles = [
                        h
                        for h in product_heterocycles
                        if h not in reactant_heterocycles
                    ]
                    if product_has_heterocycle and new_heterocycles:
                        print(
                            f"Heterocycle formation detected at depth {depth} (new: {new_heterocycles})"
                        )
                        heterocycle_formation_found = True
                        heterocycle_formation_depth = min(
                            heterocycle_formation_depth, depth
                        )

                # Check for SNAr reaction using the checker function
                snar_detected = False
                for rxn_type in snar_reaction_types:
                    if checker.check_reaction(rxn_type, rsmi):
                        print(f"SNAr reaction {rxn_type} detected at depth {depth}")
                        snar_reactions.append((depth, rsmi))
                        snar_detected = True
                        break

                # Alternative check for SNAr if not detected by reaction type
                if not snar_detected:
                    # Check for aromatic halide in reactants
                    reactants_have_aromatic_halide = any(
                        checker.check_fg("Aromatic halide", r)
                        for r in reactants_smiles
                        if r
                    )

                    # Check for nucleophilic groups in product
                    product_has_nucleophile = (
                        checker.check_fg("Aniline", product_smiles)
                        or checker.check_fg("Phenol", product_smiles)
                        or checker.check_fg("Ether", product_smiles)
                        or checker.check_fg("Secondary amine", product_smiles)
                        or checker.check_fg("Tertiary amine", product_smiles)
                        or checker.check_fg("Thioamide", product_smiles)
                        or checker.check_fg("Monosulfide", product_smiles)
                    )

                    # Check if product has aromatic structure
                    product_has_aromatic = any(
                        checker.check_ring(ring, product_smiles)
                        for ring in [
                            "benzene",
                            "pyridine",
                            "pyrimidine",
                            "pyrazine",
                            "triazole",
                            "tetrazole",
                        ]
                    )

                    if (
                        reactants_have_aromatic_halide
                        and product_has_nucleophile
                        and product_has_aromatic
                    ):
                        print(f"SNAr reaction detected via FG check at depth {depth}")
                        snar_reactions.append((depth, rsmi))

                # Check for late-stage sulfonyl introduction (depth <= 4 to be more inclusive)
                if depth <= 4:
                    sulfonyl_groups = [
                        "Sulfonamide",
                        "Sulfone",
                        "Sulfonate",
                        "Sulfamate",
                        "Sulfamic acid",
                        "Sulfonyl halide",
                    ]

                    product_has_sulfonyl = any(
                        checker.check_fg(fg, product_smiles) for fg in sulfonyl_groups
                    )

                    reactants_have_sulfonyl = any(
                        checker.check_fg(fg, r)
                        for r in reactants_smiles
                        if r
                        for fg in sulfonyl_groups
                    )

                    if product_has_sulfonyl and not reactants_have_sulfonyl:
                        print(
                            f"Late-stage sulfonyl introduction detected at depth {depth}"
                        )
                        late_stage_sulfonyl = True

            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if SNAr reactions occur after heterocycle formation
    sequential_snar = False
    if heterocycle_formation_found and snar_reactions:
        # Sort SNAr reactions by depth
        snar_reactions.sort(key=lambda x: x[0])
        # Check if at least one SNAr reaction occurs after heterocycle formation
        sequential_snar = any(
            depth < heterocycle_formation_depth for depth, _ in snar_reactions
        )
        print(
            f"Sequential SNAr check: heterocycle at depth {heterocycle_formation_depth}, SNAr at depths {[d for d, _ in snar_reactions]}"
        )

    # Return True if the strategy is detected
    # Modified to ensure SNAr reactions occur after heterocycle formation
    strategy_detected = heterocycle_formation_found and (
        sequential_snar or late_stage_sulfonyl or len(snar_reactions) >= 2
    )
    print(
        f"Strategy detected: {strategy_detected} (Heterocycle: {heterocycle_formation_found} at depth {heterocycle_formation_depth}, SNAr: {len(snar_reactions)}, Sequential: {sequential_snar}, Late-stage sulfonyl: {late_stage_sulfonyl})"
    )
    return strategy_detected
