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
    Detects if the synthesis route involves SNAr reaction for heterocycle modification.
    """
    # List of common nitrogen-containing heterocycles
    nitrogen_heterocycles = [
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
        "purine",
        "benzimidazole",
        "benzoxazole",
        "benzothiazole",
        "indazole",
        "benzotriazole",
    ]

    snar_found = False

    def dfs_traverse(node, depth=0):
        nonlocal snar_found

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            print(f"Analyzing reaction at depth {depth}: {rsmi}")

            # Check if this is a nucleophilic aromatic substitution reaction using direct reaction checks
            is_snar_direct = (
                checker.check_reaction("heteroaromatic_nuc_sub", rsmi)
                or checker.check_reaction("nucl_sub_aromatic_ortho_nitro", rsmi)
                or checker.check_reaction("nucl_sub_aromatic_para_nitro", rsmi)
            )

            if is_snar_direct:
                print(f"Direct SNAr reaction check passed: {rsmi}")
                snar_found = True
            else:
                # Check for aryl halide in reactants (leaving group)
                has_aryl_halide = any(
                    checker.check_fg("Aromatic halide", r) for r in reactants
                )

                # Check for nucleophiles in reactants
                has_nucleophile = any(
                    checker.check_fg("Primary amine", r)
                    or checker.check_fg("Secondary amine", r)
                    or checker.check_fg("Aniline", r)
                    or checker.check_fg("Primary alcohol", r)
                    or checker.check_fg("Secondary alcohol", r)
                    or checker.check_fg("Phenol", r)
                    or checker.check_fg("Aliphatic thiol", r)
                    or checker.check_fg("Aromatic thiol", r)
                    for r in reactants
                )

                # Check for nitrogen-containing heterocycle in reactants and product
                has_heterocycle_reactant = any(
                    any(checker.check_ring(ring, r) for ring in nitrogen_heterocycles)
                    for r in reactants
                )
                has_heterocycle_product = any(
                    checker.check_ring(ring, product) for ring in nitrogen_heterocycles
                )

                # Check if aromatic halide is consumed (should not be in product)
                has_aryl_halide_product = checker.check_fg("Aromatic halide", product)

                # Print diagnostic information
                print(f"  Has aryl halide in reactants: {has_aryl_halide}")
                print(f"  Has nucleophile in reactants: {has_nucleophile}")
                print(f"  Has heterocycle in reactants: {has_heterocycle_reactant}")
                print(f"  Has heterocycle in product: {has_heterocycle_product}")
                print(f"  Has aryl halide in product: {has_aryl_halide_product}")

                # SNAr pattern check: needs aryl halide and nucleophile in reactants,
                # heterocycle in either reactants or product, and aryl halide should be consumed
                is_snar_pattern = (
                    has_aryl_halide
                    and has_nucleophile
                    and (has_heterocycle_reactant or has_heterocycle_product)
                    and (
                        not has_aryl_halide_product
                        or has_aryl_halide_product != has_aryl_halide
                    )
                )

                if is_snar_pattern:
                    print(
                        f"Pattern-based SNAr heterocycle modification detected: {rsmi}"
                    )
                    snar_found = True

        # Recursively process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    print(f"SNAr heterocycle modification found: {snar_found}")
    return snar_found
