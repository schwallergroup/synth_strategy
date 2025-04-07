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
    This function detects if the synthesis concludes with an ester formation
    from a carboxylic acid (depth 0 or 1).

    Note: Since we're traversing retrosynthetically, we're actually looking for
    ester hydrolysis reactions (ester → acid) which represent ester formation
    in the forward direction.
    """
    has_late_stage_ester_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal has_late_stage_ester_formation

        # Check if this is a reaction node with metadata
        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]

            # Only check reactions at late stage (depth 0 or 1)
            if depth <= 1:
                print(f"Examining reaction at depth {depth}: {rsmi}")

                # Parse reaction components
                try:
                    reactants_part = rsmi.split(">")[0]
                    product_part = rsmi.split(">")[-1]
                    reactants = reactants_part.split(".")
                    product = product_part

                    # Check if this is an ester hydrolysis reaction (retrosynthetically represents ester formation)
                    is_ester_hydrolysis = any(
                        [
                            checker.check_reaction(
                                "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
                                rsmi,
                            ),
                            checker.check_reaction(
                                "Ester saponification (methyl deprotection)", rsmi
                            ),
                            checker.check_reaction(
                                "Ester saponification (alkyl deprotection)", rsmi
                            ),
                            checker.check_reaction("COOH ethyl deprotection", rsmi),
                        ]
                    )

                    # Also check for forward esterification reactions (in case the reaction is stored in reverse)
                    is_esterification = any(
                        [
                            checker.check_reaction("Esterification of Carboxylic Acids", rsmi),
                            checker.check_reaction("Schotten-Baumann to ester", rsmi),
                            checker.check_reaction(
                                "O-alkylation of carboxylic acids with diazo compounds", rsmi
                            ),
                            checker.check_reaction(
                                "Oxidative esterification of primary alcohols", rsmi
                            ),
                            checker.check_reaction("Transesterification", rsmi),
                            checker.check_reaction("Acetic anhydride and alcohol to ester", rsmi),
                        ]
                    )

                    print(f"  Is recognized ester hydrolysis reaction: {is_ester_hydrolysis}")
                    print(f"  Is recognized esterification reaction: {is_esterification}")

                    # If it's a recognized reaction type, we've found what we're looking for
                    if is_ester_hydrolysis or is_esterification:
                        print(f"Found late-stage ester formation at depth {depth}")
                        has_late_stage_ester_formation = True
                        return

                    # If not a recognized reaction type, check for functional group transformation
                    # In retrosynthetic direction: check if reactants contain ester and product contains acid
                    reactant_has_ester = False
                    for reactant in reactants:
                        if checker.check_fg("Ester", reactant):
                            reactant_has_ester = True
                            break

                    product_has_acid = checker.check_fg("Carboxylic acid", product)

                    # Also check the forward direction in case the reaction is stored that way
                    reactant_has_acid = False
                    for reactant in reactants:
                        if checker.check_fg("Carboxylic acid", reactant):
                            reactant_has_acid = True
                            break

                    product_has_ester = checker.check_fg("Ester", product)

                    print(f"  Reactant has ester: {reactant_has_ester}")
                    print(f"  Product has carboxylic acid: {product_has_acid}")
                    print(f"  Reactant has carboxylic acid: {reactant_has_acid}")
                    print(f"  Product has ester: {product_has_ester}")

                    # Verify the transformation by checking both directions
                    if (reactant_has_ester and product_has_acid) or (
                        reactant_has_acid and product_has_ester
                    ):
                        print(
                            f"Found late-stage ester formation (by functional group analysis) at depth {depth}"
                        )
                        has_late_stage_ester_formation = True
                        return

                except Exception as e:
                    print(f"Error analyzing reaction: {e}")

        # Recursively traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    print(f"Has late-stage ester formation: {has_late_stage_ester_formation}")
    return has_late_stage_ester_formation
