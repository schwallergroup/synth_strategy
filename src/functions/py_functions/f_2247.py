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
    This function detects a synthetic strategy that uses a late-stage amide coupling
    as the final step in the synthesis.
    """
    # Track if we find the key feature
    has_late_stage_amide_coupling = False

    def dfs_traverse(node, depth=0):
        nonlocal has_late_stage_amide_coupling

        # Check reactions at the last or second-to-last step (depth 0 or 1)
        if node["type"] == "reaction" and depth <= 1:
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                product_part = rsmi.split(">")[-1]

                reactants = reactants_part.split(".")

                print(f"Analyzing reaction at depth {depth}: {rsmi}")

                # Check if this is an amide coupling reaction
                is_amide_coupling = False

                # Check for specific amide coupling reactions
                if checker.check_reaction(
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                    rsmi,
                ):
                    is_amide_coupling = True
                    print(
                        f"Found amide coupling at depth {depth}: Acylation of Nitrogen Nucleophiles by Acyl Halides"
                    )
                elif checker.check_reaction(
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_OS",
                    rsmi,
                ):
                    is_amide_coupling = True
                    print(
                        f"Found amide coupling at depth {depth}: Acylation of Nitrogen Nucleophiles by Acyl Halides (OS)"
                    )
                elif checker.check_reaction(
                    "Acyl chloride with primary amine to amide (Schotten-Baumann)", rsmi
                ):
                    is_amide_coupling = True
                    print(
                        f"Found amide coupling at depth {depth}: Schotten-Baumann reaction"
                    )
                elif checker.check_reaction(
                    "Acyl chloride with secondary amine to amide", rsmi
                ):
                    is_amide_coupling = True
                    print(
                        f"Found amide coupling at depth {depth}: Acyl chloride with secondary amine"
                    )
                elif checker.check_reaction(
                    "Carboxylic acid with primary amine to amide", rsmi
                ):
                    is_amide_coupling = True
                    print(
                        f"Found amide coupling at depth {depth}: Carboxylic acid with amine"
                    )
                elif checker.check_reaction("Ester with primary amine to amide", rsmi):
                    is_amide_coupling = True
                    print(
                        f"Found amide coupling at depth {depth}: Ester with primary amine"
                    )
                elif checker.check_reaction(
                    "Ester with secondary amine to amide", rsmi
                ):
                    is_amide_coupling = True
                    print(
                        f"Found amide coupling at depth {depth}: Ester with secondary amine"
                    )
                elif checker.check_reaction(
                    "Acyl chloride with ammonia to amide", rsmi
                ):
                    is_amide_coupling = True
                    print(
                        f"Found amide coupling at depth {depth}: Acyl chloride with ammonia"
                    )
                elif checker.check_reaction("Ester with ammonia to amide", rsmi):
                    is_amide_coupling = True
                    print(f"Found amide coupling at depth {depth}: Ester with ammonia")
                elif checker.check_reaction("Acylation of primary amines", rsmi):
                    is_amide_coupling = True
                    print(
                        f"Found amide coupling at depth {depth}: Acylation of primary amines"
                    )
                elif checker.check_reaction("Acylation of secondary amines", rsmi):
                    is_amide_coupling = True
                    print(
                        f"Found amide coupling at depth {depth}: Acylation of secondary amines"
                    )
                elif checker.check_reaction("{Schotten-Baumann_amide}", rsmi):
                    is_amide_coupling = True
                    print(
                        f"Found amide coupling at depth {depth}: Schotten-Baumann amide"
                    )

                # If no specific reaction detected, check for reactants and products
                if not is_amide_coupling:
                    print(
                        "No specific amide coupling reaction detected, checking reactants and products..."
                    )

                    # Check for acid chloride, carboxylic acid, or ester in reactants
                    has_acyl_donor = (
                        any(checker.check_fg("Acyl halide", r) for r in reactants)
                        or any(
                            checker.check_fg("Carboxylic acid", r) for r in reactants
                        )
                        or any(checker.check_fg("Ester", r) for r in reactants)
                    )

                    # Check for amine in reactants (primary, secondary, or aniline)
                    has_amine = (
                        any(checker.check_fg("Primary amine", r) for r in reactants)
                        or any(
                            checker.check_fg("Secondary amine", r) for r in reactants
                        )
                        or any(checker.check_fg("Aniline", r) for r in reactants)
                    )

                    # Check for amide in product
                    has_amide_product = (
                        checker.check_fg("Primary amide", product_part)
                        or checker.check_fg("Secondary amide", product_part)
                        or checker.check_fg("Tertiary amide", product_part)
                    )

                    print(f"Acyl donor present: {has_acyl_donor}")
                    print(f"Amine present: {has_amine}")
                    print(f"Amide in product: {has_amide_product}")

                    # Check if we have the necessary components for amide coupling
                    if has_acyl_donor and has_amine and has_amide_product:
                        # Additional check: verify that a new amide bond is formed
                        # This is a simplified check - in a real implementation, we would need to
                        # track atom mappings to ensure the specific atoms form a new bond
                        is_amide_coupling = True
                        print(
                            f"Found amide coupling at depth {depth}: Generic amide formation"
                        )

                if is_amide_coupling:
                    has_late_stage_amide_coupling = True
            except Exception as e:
                print(f"Error analyzing reaction at depth {depth}: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    if has_late_stage_amide_coupling:
        print("Detected strategy: Late-stage amide coupling")
    else:
        print("Strategy not detected: No late-stage amide coupling found")

    return has_late_stage_amide_coupling
