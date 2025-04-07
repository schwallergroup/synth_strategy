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
    Detects if the final step (depth 1) is an amide coupling reaction
    """
    final_step_is_amide_coupling = False

    def dfs_traverse(node, depth=0):
        nonlocal final_step_is_amide_coupling

        print(
            f"Traversing node at depth {depth}, type: {node['type']}, smiles: {node.get('smiles', 'N/A')}"
        )

        # Check if this is a reaction node at the final step (depth 1)
        if node["type"] == "reaction" and depth == 1:
            try:
                rsmi = node["metadata"]["rsmi"]
                print(
                    f"Analyzing potential final step reaction at depth {depth}: {rsmi}"
                )

                # Check if this is an amide coupling using reaction checkers
                amide_coupling_reactions = [
                    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                    "Carboxylic acid with primary amine to amide",
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                    "Acyl chloride with secondary amine to amide",
                    "Ester with primary amine to amide",
                    "Ester with secondary amine to amide",
                    "Schotten-Baumann_amide",
                    "Acylation of primary amines",
                    "Acylation of secondary amines",
                    "Acyl chloride with ammonia to amide",
                    "Ester with ammonia to amide",
                ]

                for reaction_type in amide_coupling_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        print(f"Found amide coupling as final step: {reaction_type}")
                        final_step_is_amide_coupling = True
                        return

                # If reaction check fails, try checking for functional group changes
                if not final_step_is_amide_coupling:
                    print("Reaction check failed, trying functional group analysis")

                    reactants_part = rsmi.split(">")[0]
                    product_part = rsmi.split(">")[-1]
                    reactants = reactants_part.split(".")

                    # Check for carboxylic acid derivatives in reactants
                    has_carboxylic_acid = any(
                        checker.check_fg("Carboxylic acid", r) for r in reactants
                    )
                    has_acyl_halide = any(
                        checker.check_fg("Acyl halide", r) for r in reactants
                    )
                    has_ester = any(checker.check_fg("Ester", r) for r in reactants)
                    has_anhydride = any(
                        checker.check_fg("Anhydride", r) for r in reactants
                    )

                    print(
                        f"Carboxylic acid: {has_carboxylic_acid}, Acyl halide: {has_acyl_halide}, Ester: {has_ester}, Anhydride: {has_anhydride}"
                    )

                    # Check for amine in reactants
                    has_primary_amine = any(
                        checker.check_fg("Primary amine", r) for r in reactants
                    )
                    has_secondary_amine = any(
                        checker.check_fg("Secondary amine", r) for r in reactants
                    )
                    has_aniline = any(checker.check_fg("Aniline", r) for r in reactants)

                    print(
                        f"Primary amine: {has_primary_amine}, Secondary amine: {has_secondary_amine}, Aniline: {has_aniline}"
                    )

                    # Check for amide in product
                    has_primary_amide = checker.check_fg("Primary amide", product_part)
                    has_secondary_amide = checker.check_fg(
                        "Secondary amide", product_part
                    )
                    has_tertiary_amide = checker.check_fg(
                        "Tertiary amide", product_part
                    )

                    print(
                        f"Product - Primary amide: {has_primary_amide}, Secondary amide: {has_secondary_amide}, Tertiary amide: {has_tertiary_amide}"
                    )

                    # Check if we have the necessary components for amide coupling
                    has_acid_derivative = (
                        has_carboxylic_acid
                        or has_acyl_halide
                        or has_ester
                        or has_anhydride
                    )
                    has_amine = has_primary_amine or has_secondary_amine or has_aniline
                    has_amide = (
                        has_primary_amide or has_secondary_amide or has_tertiary_amide
                    )

                    if has_acid_derivative and has_amine and has_amide:
                        print(
                            f"Found amide coupling as final step (via functional group analysis)"
                        )
                        final_step_is_amide_coupling = True
                        return
            except Exception as e:
                print(f"Error analyzing reaction: {e}")

        # Traverse children with incremented depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    print(f"Final result: {final_step_is_amide_coupling}")

    return final_step_is_amide_coupling
