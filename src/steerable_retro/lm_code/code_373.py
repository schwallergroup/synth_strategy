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
    This function detects a synthetic strategy involving late-stage amide formation
    (amide formation in the final step of the synthesis).
    """
    # Define amide formation reaction types
    amide_formation_reactions = [
        "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
        "Carboxylic acid with primary amine to amide",
        "Ester with primary amine to amide",
        "Ester with secondary amine to amide",
        "Ester with ammonia to amide",
        "Acyl chloride with primary amine to amide (Schotten-Baumann)",
        "Acyl chloride with secondary amine to amide",
        "Acyl chloride with ammonia to amide",
        "Schotten-Baumann_amide",
        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
        "Aminolysis of esters",
    ]

    # Track if we found amide formation at the first reaction (late stage)
    found_late_stage_amide = [False]

    def traverse_route(node, depth=0):
        # If we already found a late-stage amide formation, no need to continue
        if found_late_stage_amide[0]:
            return

        # For reaction nodes, check if it's an amide formation reaction
        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            print(f"Checking reaction at depth {depth}: {rsmi}")

            # Consider both depth 0 and depth 1 as late-stage reactions
            if depth <= 1:
                print(f"Evaluating late-stage reaction at depth {depth}")
                # Check for amide formation using the checker function
                for rxn_name in amide_formation_reactions:
                    if checker.check_reaction(rxn_name, rsmi):
                        print(f"Found late-stage amide formation: {rxn_name}")

                        # Verify by checking product for amide group
                        product = rsmi.split(">")[-1]
                        if (
                            checker.check_fg("Primary amide", product)
                            or checker.check_fg("Secondary amide", product)
                            or checker.check_fg("Tertiary amide", product)
                        ):
                            found_late_stage_amide[0] = True
                            print(f"Confirmed amide group in product: {product}")
                            return

                # If no direct reaction match, check if the product contains an amide
                # and reactants don't (indicating amide formation)
                product = rsmi.split(">")[-1]
                reactants = rsmi.split(">")[0].split(".")

                product_has_amide = (
                    checker.check_fg("Primary amide", product)
                    or checker.check_fg("Secondary amide", product)
                    or checker.check_fg("Tertiary amide", product)
                )

                reactants_have_amide = any(
                    checker.check_fg("Primary amide", reactant)
                    or checker.check_fg("Secondary amide", reactant)
                    or checker.check_fg("Tertiary amide", reactant)
                    for reactant in reactants
                )

                if product_has_amide and not reactants_have_amide:
                    print(f"Detected amide formation by functional group analysis")
                    found_late_stage_amide[0] = True
                    return

        # Recursively check children with increased depth
        for child in node.get("children", []):
            traverse_route(child, depth + 1)

    # Start traversal from the root
    traverse_route(route)

    print(f"Late-stage amide formation detection result: {found_late_stage_amide[0]}")
    return found_late_stage_amide[0]
