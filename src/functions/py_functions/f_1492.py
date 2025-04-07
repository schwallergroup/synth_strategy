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
    This function detects aromatic sulfonation (introduction of sulfonyl group to aromatic ring).
    Aromatic sulfonation typically involves the introduction of a sulfonic acid (SO3H) group
    to an aromatic ring, or related transformations involving sulfonyl groups.
    """
    found_sulfonation = False

    def dfs_traverse(node):
        nonlocal found_sulfonation

        if found_sulfonation:
            return

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for aromatic sulfonation using the checker function
            if checker.check_reaction("Aromatic sulfonyl chlorination", rsmi):
                print(f"Found aromatic sulfonation via reaction check: {rsmi}")
                found_sulfonation = True
                return

            # Check for sulfonic acid or sulfonate formation on aromatic rings
            has_aromatic_in_reactants = any(
                checker.check_ring("benzene", reactant) for reactant in reactants
            )
            has_aromatic_in_product = checker.check_ring("benzene", product)

            # Check for sulfonic acid or sulfonate in product
            has_sulfonic_acid_in_product = checker.check_fg(
                "Sulfonic acid", product
            ) or checker.check_fg("Sulfonate", product)

            # Check for sulfonyl halide in reactants (could be used for sulfonation)
            has_sulfonyl_halide_in_reactants = any(
                checker.check_fg("Sulfonyl halide", reactant) for reactant in reactants
            )

            # Check for sulfonamide formation (a related sulfonation process)
            has_sulfonamide_in_product = checker.check_fg("Sulfonamide", product)

            # Check for sulfonate esters
            has_sulfonate_in_product = checker.check_fg("Sulfonate", product)

            # Check for SO3 or related reagents in reactants
            has_so3_reagent = any(
                "SO3" in reactant or "SO2" in reactant for reactant in reactants
            )

            # Various sulfonation scenarios
            if has_aromatic_in_reactants and has_aromatic_in_product:
                if has_sulfonic_acid_in_product and not any(
                    checker.check_fg("Sulfonic acid", reactant)
                    for reactant in reactants
                ):
                    print(
                        f"Found aromatic sulfonation (sulfonic acid formation): {rsmi}"
                    )
                    found_sulfonation = True
                    return

                if has_sulfonate_in_product and not any(
                    checker.check_fg("Sulfonate", reactant) for reactant in reactants
                ):
                    print(f"Found aromatic sulfonation (sulfonate formation): {rsmi}")
                    found_sulfonation = True
                    return

                if has_sulfonyl_halide_in_reactants and has_sulfonamide_in_product:
                    print(f"Found aromatic sulfonamide formation: {rsmi}")
                    found_sulfonation = True
                    return

                if has_so3_reagent and (
                    has_sulfonic_acid_in_product or has_sulfonate_in_product
                ):
                    print(f"Found direct aromatic sulfonation with SO3: {rsmi}")
                    found_sulfonation = True
                    return

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return found_sulfonation
