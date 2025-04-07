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
    Detects if the synthesis follows a linear strategy with a late-stage C-C bond formation
    (typically via Grignard or similar reaction).
    """
    # Track if we found a late-stage C-C bond formation
    late_stage_cc_bond = False

    # Track reaction steps and their depths to determine linearity
    reaction_depths = []

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_cc_bond

        if node["type"] == "reaction":
            # Store the depth of this reaction
            reaction_depths.append(depth)

            # Check if this is a late-stage reaction (depth 0 or 1)
            if depth <= 1 and "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]

                # Check for C-C bond forming reactions
                cc_bond_forming_reactions = [
                    "Grignard from aldehyde to alcohol",
                    "Grignard from ketone to alcohol",
                    "Suzuki coupling with boronic acids",
                    "Suzuki coupling with boronic esters",
                    "Suzuki coupling with boronic acids OTf",
                    "Suzuki coupling with boronic esters OTf",
                    "Suzuki coupling with sulfonic esters",
                    "Negishi coupling",
                    "Wittig reaction with triphenylphosphorane",
                    "Wittig with Phosphonium",
                    "Heck terminal vinyl",
                    "Oxidative Heck reaction",
                    "Oxidative Heck reaction with vinyl ester",
                    "Stille reaction_aryl",
                    "Stille reaction_vinyl",
                    "Stille reaction_benzyl",
                    "Stille reaction_allyl",
                    "Kumada cross-coupling",
                    "Aryllithium cross-coupling",
                    "Friedel-Crafts alkylation",
                    "Friedel-Crafts alkylation with halide",
                    "Knoevenagel Condensation",
                    "Aldol condensation",
                    "Michael addition",
                    "Michael addition methyl",
                    "Hiyama-Denmark Coupling",
                    "decarboxylative_coupling",
                    "Catellani reaction ortho",
                    "Catellani reaction para",
                    "beta C(sp3) arylation",
                ]

                for reaction_type in cc_bond_forming_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        print(
                            f"Found late-stage C-C bond formation via {reaction_type}"
                        )
                        late_stage_cc_bond = True
                        break

                # If no specific reaction type matched, check for general C-C bond formation
                if not late_stage_cc_bond:
                    try:
                        # Extract reactants and product
                        reactants_str = rsmi.split(">")[0]
                        product_str = rsmi.split(">")[-1]

                        # Check if this is a Grignard reaction not caught by the specific types
                        for reactant in reactants_str.split("."):
                            if checker.check_fg("Magnesium halide", reactant):
                                print(
                                    "Found late-stage C-C bond formation with Grignard reagent"
                                )
                                late_stage_cc_bond = True
                                break

                            # Check for organolithium compounds
                            if checker.check_fg(
                                "Alkyl lithium", reactant
                            ) or checker.check_fg("Aryl lithium", reactant):
                                print(
                                    "Found late-stage C-C bond formation with organolithium reagent"
                                )
                                late_stage_cc_bond = True
                                break
                    except Exception as e:
                        print(f"Error analyzing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if the synthesis is linear
    is_linear = len(reaction_depths) >= 2  # At least 2 steps

    if reaction_depths:
        # Sort depths to analyze the pattern
        reaction_depths.sort()

        # Check if the depths follow a pattern consistent with linear synthesis
        # For a linear synthesis, depths should increase with a consistent pattern
        # Allow gaps of up to 2 between consecutive depths
        max_gap = 0
        for i in range(1, len(reaction_depths)):
            gap = reaction_depths[i] - reaction_depths[i - 1]
            max_gap = max(max_gap, gap)

        # If max gap between consecutive depths is > 2, it suggests significant branching
        is_linear = is_linear and max_gap <= 2

        # Additional check: in a linear synthesis, the number of unique depths should be close to the total number of reactions
        unique_depths = len(set(reaction_depths))
        is_linear = (
            is_linear and unique_depths >= len(reaction_depths) * 0.7
        )  # At least 70% of reactions should be at unique depths

    print(f"Reaction depths: {reaction_depths}")
    print(f"Is linear: {is_linear}, Has late-stage C-C bond: {late_stage_cc_bond}")

    # Return True if it's a linear synthesis with late-stage C-C bond formation
    return is_linear and late_stage_cc_bond
