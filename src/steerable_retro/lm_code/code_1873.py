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
    Detects a synthetic strategy involving benzylic halide activation of a heterocyclic core
    followed by nucleophilic substitution with a nitrogen-containing fragment in the late stage
    of the synthesis.
    """
    benzylic_halide_found = False
    late_stage_coupling = False

    print("Starting analysis for benzylic halide coupling strategy...")

    def dfs_traverse(node, depth=0):
        nonlocal benzylic_halide_found, late_stage_coupling

        if node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"\nAnalyzing reaction at depth {depth}: {rsmi}")

                # Check if this is a nucleophilic substitution reaction
                is_nucleophilic_sub = False

                # Check for specific nucleophilic substitution reactions
                nucleophilic_reactions = [
                    "N-alkylation of primary amines with alkyl halides",
                    "N-alkylation of secondary amines with alkyl halides",
                    "Williamson Ether Synthesis",
                    "S-alkylation of thiols",
                    "Buchwald-Hartwig",
                    "N-arylation",
                    "Mitsunobu",
                    "SN2",
                ]

                for rxn_type in nucleophilic_reactions:
                    if checker.check_reaction(rxn_type, rsmi):
                        is_nucleophilic_sub = True
                        print(f"Found nucleophilic substitution reaction: {rxn_type}")
                        break

                # If no specific reaction type matched, check for general substitution pattern
                if not is_nucleophilic_sub:
                    # Look for a pattern where a halide is replaced by a nitrogen
                    halide_reactant = None
                    nitrogen_reactant = None

                    for reactant in reactants:
                        if (
                            checker.check_fg("Primary halide", reactant)
                            or checker.check_fg("Secondary halide", reactant)
                            or checker.check_fg("Tertiary halide", reactant)
                        ):
                            halide_reactant = reactant

                        if (
                            checker.check_fg("Primary amine", reactant)
                            or checker.check_fg("Secondary amine", reactant)
                            or checker.check_fg("Tertiary amine", reactant)
                            or checker.check_fg("Aniline", reactant)
                        ):
                            nitrogen_reactant = reactant

                    if halide_reactant and nitrogen_reactant:
                        # Check if product has new C-N bond and no halide
                        has_halide_in_product = (
                            checker.check_fg("Primary halide", product)
                            or checker.check_fg("Secondary halide", product)
                            or checker.check_fg("Tertiary halide", product)
                        )

                        if not has_halide_in_product:
                            is_nucleophilic_sub = True
                            print("Found nucleophilic substitution pattern (halide replacement)")

                if is_nucleophilic_sub:
                    print(f"Found nucleophilic substitution reaction at depth {depth}")

                    # Identify benzylic halide reactant and nitrogen nucleophile
                    benzylic_halide_reactant = None
                    nitrogen_nucleophile_reactant = None

                    for reactant in reactants:
                        # Check for halide
                        has_halide = (
                            checker.check_fg("Primary halide", reactant)
                            or checker.check_fg("Secondary halide", reactant)
                            or checker.check_fg("Tertiary halide", reactant)
                        )

                        if has_halide:
                            # Check if it's benzylic (connected to aromatic carbon)
                            mol = Chem.MolFromSmiles(reactant)
                            if mol:
                                # More comprehensive benzylic halide check
                                benzylic_pattern = Chem.MolFromSmarts("c[CD1,CD2,CD3][F,Cl,Br,I]")
                                if mol.HasSubstructMatch(benzylic_pattern):
                                    print(f"Found potential benzylic halide: {reactant}")

                                    # Check if the core is heterocyclic
                                    heterocycle_found = False
                                    heterocycle_name = None

                                    # Check for common heterocycles
                                    heterocycles = [
                                        "pyridine",
                                        "pyrimidine",
                                        "pyrazine",
                                        "pyridazine",
                                        "imidazole",
                                        "thiazole",
                                        "oxazole",
                                        "triazole",
                                        "tetrazole",
                                        "indole",
                                        "benzimidazole",
                                        "quinoline",
                                        "isoquinoline",
                                        "furan",
                                        "thiophene",
                                        "pyrrole",
                                        "oxadiazole",
                                        "thiadiazole",
                                        "piperidine",
                                        "piperazine",
                                        "morpholine",
                                        "thiomorpholine",
                                    ]

                                    for ring in heterocycles:
                                        if checker.check_ring(ring, reactant):
                                            heterocycle_found = True
                                            heterocycle_name = ring
                                            print(f"Found heterocycle ({ring}) in benzylic halide")
                                            break

                                    if heterocycle_found:
                                        benzylic_halide_reactant = reactant

                        # Check for nitrogen nucleophile
                        has_amine = (
                            checker.check_fg("Primary amine", reactant)
                            or checker.check_fg("Secondary amine", reactant)
                            or checker.check_fg("Tertiary amine", reactant)
                            or checker.check_fg("Aniline", reactant)
                        )

                        if has_amine:
                            print(f"Found nitrogen nucleophile: {reactant}")
                            nitrogen_nucleophile_reactant = reactant

                    # Verify both components are found and the reaction is at a late stage
                    if benzylic_halide_reactant and nitrogen_nucleophile_reactant:
                        print(
                            f"Found both benzylic halide and nitrogen nucleophile at depth {depth}"
                        )

                        # Verify the product contains the heterocycle but not the halide
                        product_mol = Chem.MolFromSmiles(product)
                        if product_mol:
                            # Check if product still has the heterocycle
                            heterocycle_in_product = False
                            for ring in heterocycles:
                                if checker.check_ring(ring, product):
                                    heterocycle_in_product = True
                                    print(f"Found heterocycle ({ring}) in product")
                                    break

                            # Check if halide is gone in product
                            halide_in_product = (
                                checker.check_fg("Primary halide", product)
                                or checker.check_fg("Secondary halide", product)
                                or checker.check_fg("Tertiary halide", product)
                            )

                            if halide_in_product:
                                print("Halide still present in product")

                            if heterocycle_in_product and not halide_in_product:
                                benzylic_halide_found = True
                                print("Confirmed benzylic halide coupling transformation")

                                # Check if this is a late stage reaction (depth <= 3)
                                if depth <= 3:
                                    late_stage_coupling = True
                                    print(
                                        f"SUCCESS: Found late-stage benzylic halide coupling at depth {depth}"
                                    )
                                else:
                                    print(
                                        f"Found benzylic halide coupling, but not at late stage (depth {depth} > 3)"
                                    )
                            else:
                                print("Product structure doesn't match expected transformation")
            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(
        f"\nFinal result: benzylic_halide_found={benzylic_halide_found}, late_stage_coupling={late_stage_coupling}"
    )
    return benzylic_halide_found and late_stage_coupling
