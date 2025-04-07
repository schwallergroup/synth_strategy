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
    Detects a linear synthesis strategy with amide formation from carboxylic acid.
    """
    has_amide_formation = False
    reaction_count = 0
    max_branch_factor = 0

    def dfs_traverse(node, depth=0):
        nonlocal has_amide_formation, reaction_count, max_branch_factor

        if node["type"] == "reaction":
            reaction_count += 1

            # Check branching factor at this reaction node
            if "children" in node:
                branch_factor = len(node["children"])
                max_branch_factor = max(max_branch_factor, branch_factor)

            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Analyzing reaction at depth {depth}: {rsmi}")

                # Check for amide formation using reaction checkers
                amide_formation_reactions = [
                    "Carboxylic acid with primary amine to amide",
                    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                    "Acylation of primary amines",
                    "Acylation of secondary amines",
                    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                    "Acyl chloride with secondary amine to amide",
                    "Acyl chloride with ammonia to amide",
                    "Ester with primary amine to amide",
                    "Ester with secondary amine to amide",
                    "Ester with ammonia to amide",
                    "Schotten-Baumann_amide",
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                    "Aminolysis of esters",
                    "Acylation of secondary amines with anhydrides",
                ]

                # Check if any of the amide formation reactions match
                for rxn in amide_formation_reactions:
                    if checker.check_reaction(rxn, rsmi):
                        print(f"Matched reaction type: {rxn}")

                        # Verify the reactants have the expected functional groups
                        acid_groups = ["Carboxylic acid", "Acyl halide", "Ester", "Anhydride"]
                        amine_groups = [
                            "Primary amine",
                            "Secondary amine",
                            "Tertiary amine",
                            "Aniline",
                        ]
                        amide_groups = ["Primary amide", "Secondary amide", "Tertiary amide"]

                        has_acid_precursor = any(
                            any(checker.check_fg(fg, r) for fg in acid_groups) for r in reactants
                        )

                        has_amine = any(
                            any(checker.check_fg(fg, r) for fg in amine_groups) for r in reactants
                        )

                        has_amide_product = any(
                            checker.check_fg(fg, product) for fg in amide_groups
                        )

                        print(f"Has acid precursor: {has_acid_precursor}")
                        print(f"Has amine: {has_amine}")
                        print(f"Has amide product: {has_amide_product}")

                        # Check if reactants don't already have amide groups
                        reactants_have_amide = any(
                            any(checker.check_fg(fg, r) for fg in amide_groups) for r in reactants
                        )

                        # If reactants have necessary groups and product has amide
                        if (
                            has_acid_precursor
                            and has_amine
                            and has_amide_product
                            and not reactants_have_amide
                        ):
                            has_amide_formation = True
                            print(f"Confirmed amide formation at depth {depth}")
                            break

                # If no specific reaction type matched, try a more general approach
                if not has_amide_formation:
                    # Check if reactants have acid/amine and product has amide (that wasn't in reactants)
                    acid_groups = ["Carboxylic acid", "Acyl halide", "Ester", "Anhydride"]
                    amine_groups = ["Primary amine", "Secondary amine", "Tertiary amine", "Aniline"]
                    amide_groups = ["Primary amide", "Secondary amide", "Tertiary amide"]

                    has_acid_precursor = any(
                        any(checker.check_fg(fg, r) for fg in acid_groups) for r in reactants
                    )

                    has_amine = any(
                        any(checker.check_fg(fg, r) for fg in amine_groups) for r in reactants
                    )

                    has_amide_product = any(checker.check_fg(fg, product) for fg in amide_groups)

                    # Check if reactants don't already have amide groups
                    reactants_have_amide = any(
                        any(checker.check_fg(fg, r) for fg in amide_groups) for r in reactants
                    )

                    if (
                        has_acid_precursor
                        and has_amine
                        and has_amide_product
                        and not reactants_have_amide
                    ):
                        print("Detected amide formation through functional group analysis")
                        has_amide_formation = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if it's a linear synthesis (minimal branching and at least 3 reactions)
    is_linear = reaction_count >= 3 and max_branch_factor <= 2

    print(f"Reaction count: {reaction_count}, Max branch factor: {max_branch_factor}")
    print(f"Has amide formation: {has_amide_formation}, Is linear: {is_linear}")

    return has_amide_formation and is_linear
