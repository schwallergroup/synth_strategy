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
    Detects if the synthetic route contains amide formation in the late stage (depth 0-1).
    """
    late_stage_amide_found = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_amide_found

        if node["type"] == "reaction" and depth <= 1:
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Checking reaction at depth {depth}: {rsmi}")

                # Check for amide formation reaction types
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
                    "Carboxylic acid to amide conversion",
                ]

                # First check if it's directly an amide formation reaction
                is_amide_reaction = False
                for reaction_type in amide_formation_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        print(f"Found amide formation reaction: {reaction_type}")
                        is_amide_reaction = True
                        break

                # If not a known reaction type, check for functional group patterns
                if not is_amide_reaction:
                    # Check for amide in product
                    has_amide_product = any(
                        checker.check_fg(amide_type, product)
                        for amide_type in ["Primary amide", "Secondary amide", "Tertiary amide"]
                    )

                    if has_amide_product:
                        print(f"Found amide in product: {product}")

                        # Check for acid derivatives and amines in reactants
                        has_acid_derivative = False
                        has_amine = False

                        acid_derivatives = ["Carboxylic acid", "Acyl halide", "Ester", "Anhydride"]
                        amine_types = ["Primary amine", "Secondary amine", "Aniline"]

                        for reactant in reactants:
                            if any(checker.check_fg(acid, reactant) for acid in acid_derivatives):
                                print(f"Found acid derivative in reactant: {reactant}")
                                has_acid_derivative = True

                            if any(checker.check_fg(amine, reactant) for amine in amine_types):
                                print(f"Found amine in reactant: {reactant}")
                                has_amine = True

                        if has_acid_derivative and has_amine:
                            print(
                                f"Detected amide formation based on functional groups at depth {depth}"
                            )
                            is_amide_reaction = True

                # If we've confirmed this is an amide formation reaction
                if is_amide_reaction:
                    # Double-check that the product actually contains an amide
                    if any(
                        checker.check_fg(amide_type, product)
                        for amide_type in ["Primary amide", "Secondary amide", "Tertiary amide"]
                    ):
                        print(f"Confirmed late-stage amide formation at depth {depth}")
                        late_stage_amide_found = True
                        return  # Found what we're looking for, no need to continue

        # Traverse children with increased depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    print("Starting traversal to find late-stage amide formation...")
    dfs_traverse(route)

    if not late_stage_amide_found:
        print("No late-stage amide formation found in the route")

    return late_stage_amide_found
