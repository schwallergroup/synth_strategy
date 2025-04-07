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
    Detects if the synthesis route uses a late-stage amide formation strategy,
    where an amide bond is formed in the final step of the synthesis.
    """
    final_amide_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal final_amide_formation

        print(f"Traversing node at depth {depth}, type: {node['type']}")

        # Check if this is a reaction node
        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            reactants = reactants_part.split(".")
            product = rsmi.split(">")[-1]

            print(f"Found reaction at depth {depth}: {rsmi}")

            # Check if this is the first reaction (depth 1 in retrosynthetic route)
            # The first reaction is at depth 1 because depth 0 is the final product
            if depth == 1:
                print(f"Checking final step reaction: {rsmi}")

                # First check if this is a known amide formation reaction
                amide_formation_reactions = [
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                    "Acyl chloride with ammonia to amide",
                    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                    "Acyl chloride with secondary amine to amide",
                    "Carboxylic acid with primary amine to amide",
                    "Ester with ammonia to amide",
                    "Ester with primary amine to amide",
                    "Ester with secondary amine to amide",
                    "Schotten-Baumann_amide",
                    "Acylation of primary amines",
                    "Acylation of secondary amines",
                    "Carboxylic acid to amide conversion",
                ]

                for reaction_type in amide_formation_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        print(f"Detected amide formation reaction: {reaction_type}")
                        final_amide_formation = True
                        return

                # If no specific reaction type matched, check for amide formation by functional group analysis

                # Check if product contains an amide group
                product_has_primary_amide = checker.check_fg("Primary amide", product)
                product_has_secondary_amide = checker.check_fg(
                    "Secondary amide", product
                )
                product_has_tertiary_amide = checker.check_fg("Tertiary amide", product)

                product_has_amide = (
                    product_has_primary_amide
                    or product_has_secondary_amide
                    or product_has_tertiary_amide
                )

                if product_has_amide:
                    print("Product contains amide group")

                    # Check if at least one reactant does NOT have the same amide
                    # This confirms the amide was formed, not just preserved
                    reactants_all_have_same_amide = True

                    if product_has_primary_amide:
                        reactants_all_have_same_amide = all(
                            checker.check_fg("Primary amide", r) for r in reactants
                        )
                    elif product_has_secondary_amide:
                        reactants_all_have_same_amide = all(
                            checker.check_fg("Secondary amide", r) for r in reactants
                        )
                    elif product_has_tertiary_amide:
                        reactants_all_have_same_amide = all(
                            checker.check_fg("Tertiary amide", r) for r in reactants
                        )

                    if not reactants_all_have_same_amide:
                        print(
                            "Amide was formed in this reaction (not present in all reactants)"
                        )

                        # Check for reactants that can form amides
                        has_acid = any(
                            checker.check_fg("Carboxylic acid", reactant)
                            for reactant in reactants
                        )
                        has_acyl_halide = any(
                            checker.check_fg("Acyl halide", reactant)
                            for reactant in reactants
                        )
                        has_ester = any(
                            checker.check_fg("Ester", reactant)
                            for reactant in reactants
                        )
                        has_anhydride = any(
                            checker.check_fg("Anhydride", reactant)
                            for reactant in reactants
                        )

                        has_primary_amine = any(
                            checker.check_fg("Primary amine", reactant)
                            for reactant in reactants
                        )
                        has_secondary_amine = any(
                            checker.check_fg("Secondary amine", reactant)
                            for reactant in reactants
                        )
                        has_ammonia = "N" in reactants_part and not any(
                            atom
                            for atom in "BCDEFGHIJKLOPQRSTUVWXYZ"
                            if atom in reactants_part
                        )

                        # Check if we have appropriate reactants for amide formation
                        if (
                            has_acid or has_acyl_halide or has_ester or has_anhydride
                        ) and (has_primary_amine or has_secondary_amine or has_ammonia):
                            print("Detected appropriate reactants for amide formation")
                            final_amide_formation = True
                            return

        # Continue traversal with increased depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    print("Starting traversal of synthesis route")
    dfs_traverse(route)
    print(f"Final result: {final_amide_formation}")
    return final_amide_formation
