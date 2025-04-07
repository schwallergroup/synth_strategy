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
    Detects a synthetic strategy involving benzoxazole formation and late-stage nitro reduction.
    The strategy involves:
    1. Benzoxazole ring formation in an intermediate step
    2. Nitro reduction as the final step
    """
    benzoxazole_formed = False
    late_nitro_reduction = False

    def dfs_traverse(node, depth=0):
        nonlocal benzoxazole_formed, late_nitro_reduction

        if node["type"] == "reaction":
            # Extract reactants and product from reaction SMILES
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if not rsmi:
                return

            # Check for benzoxazole formation
            # Check if any of the benzoxazole formation reactions are detected
            benzoxazole_rxn_types = [
                "Benzoxazole formation from aldehyde",
                "Benzoxazole formation from acyl halide",
                "Benzoxazole formation from ester/carboxylic acid",
                "Benzoxazole formation (intramolecular)",
                "{benzoxazole_arom-aldehyde}",
                "{benzoxazole_carboxylic-acid}",
            ]

            for rxn_type in benzoxazole_rxn_types:
                if checker.check_reaction(rxn_type, rsmi):
                    print(f"Benzoxazole formation detected at depth: {depth}")
                    benzoxazole_formed = True
                    break

            # If reaction type check failed, try structural check
            if not benzoxazole_formed:
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                try:
                    product_has_benzoxazole = checker.check_ring("benzoxazole", product_smiles)
                    reactant_has_benzoxazole = any(
                        checker.check_ring("benzoxazole", r) for r in reactants_smiles
                    )

                    if product_has_benzoxazole and not reactant_has_benzoxazole:
                        print(f"Benzoxazole formation detected at depth: {depth}")
                        benzoxazole_formed = True
                except Exception as e:
                    print(f"Error checking benzoxazole structure: {e}")

            # Check for nitro reduction at depth 0 or 1 (final or penultimate step)
            if depth <= 1:
                # Check for nitro reduction reaction
                if checker.check_reaction("Reduction of nitro groups to amines", rsmi):
                    print(f"Late-stage nitro reduction detected at depth: {depth}")
                    late_nitro_reduction = True
                else:
                    # Fallback to structural check
                    reactants_smiles = rsmi.split(">")[0].split(".")
                    product_smiles = rsmi.split(">")[-1]

                    try:
                        reactant_has_nitro = any(
                            checker.check_fg("Nitro group", r) for r in reactants_smiles
                        )
                        # Check for both primary amine and aniline in product
                        product_has_amine = checker.check_fg(
                            "Primary amine", product_smiles
                        ) or checker.check_fg("Aniline", product_smiles)
                        # Ensure nitro group is removed in product
                        product_has_no_nitro = not checker.check_fg("Nitro group", product_smiles)

                        if reactant_has_nitro and product_has_amine and product_has_no_nitro:
                            print(f"Late-stage nitro reduction detected at depth: {depth}")
                            late_nitro_reduction = True
                    except Exception as e:
                        print(f"Error checking nitro reduction: {e}")

        # Process children with incremented depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Return True if both conditions are met
    return benzoxazole_formed and late_nitro_reduction
