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
    Detects if the synthesis route involves a late-stage amide coupling
    (formation of an amide bond in the final or penultimate step).
    """
    found_late_amide = False
    max_depth = 1  # Consider only reactions at depth 0 or 1 as "late-stage"

    def dfs_traverse(node, current_depth=0):
        nonlocal found_late_amide

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]

                # Only check reactions at depth 0 or 1 (late-stage)
                if current_depth <= max_depth:
                    print(f"Checking reaction at depth {current_depth}: {rsmi}")

                    # Check for various amide coupling reactions
                    amide_reactions = [
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
                    ]

                    for reaction_type in amide_reactions:
                        if checker.check_reaction(reaction_type, rsmi):
                            print(
                                f"Found late-stage amide coupling: {reaction_type} at depth {current_depth}"
                            )
                            found_late_amide = True
                            break

                    # Additional check for amide formation if no specific reaction type matched
                    if not found_late_amide:
                        try:
                            reactants = rsmi.split(">")[0].split(".")
                            product = rsmi.split(">")[-1]

                            # Check if any reactant has carboxylic acid or acyl halide
                            has_carboxylic_acid = any(
                                checker.check_fg("Carboxylic acid", r)
                                for r in reactants
                            )
                            has_acyl_halide = any(
                                checker.check_fg("Acyl halide", r) for r in reactants
                            )

                            # Check if any reactant has amine
                            has_primary_amine = any(
                                checker.check_fg("Primary amine", r) for r in reactants
                            )
                            has_secondary_amine = any(
                                checker.check_fg("Secondary amine", r)
                                for r in reactants
                            )

                            # Check if product has amide
                            has_primary_amide = checker.check_fg(
                                "Primary amide", product
                            )
                            has_secondary_amide = checker.check_fg(
                                "Secondary amide", product
                            )
                            has_tertiary_amide = checker.check_fg(
                                "Tertiary amide", product
                            )

                            # Verify amide formation
                            if (
                                (has_carboxylic_acid or has_acyl_halide)
                                and (has_primary_amine or has_secondary_amine)
                                and (
                                    has_primary_amide
                                    or has_secondary_amide
                                    or has_tertiary_amide
                                )
                            ):
                                print(
                                    f"Found late-stage amide coupling (FG analysis) at depth {current_depth}"
                                )
                                found_late_amide = True
                        except Exception as e:
                            print(f"Error analyzing reaction: {e}")

        # Traverse children with incremented depth
        for child in node.get("children", []):
            if (
                not found_late_amide
            ):  # Stop traversal if we already found what we're looking for
                dfs_traverse(child, current_depth + 1)

    # Start traversal
    dfs_traverse(route)

    return found_late_amide
