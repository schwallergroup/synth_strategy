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
    This function detects if the final step (depth 1) in the synthesis is an amide formation.
    """
    final_step_is_amide_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal final_step_is_amide_formation

        print(
            f"Traversing node at depth {depth}, type: {node['type']}, smiles: {node.get('smiles', 'N/A')}"
        )

        # Check if this is a reaction node at the final step (depth 1)
        if node["type"] == "reaction" and depth == 1:
            try:
                rsmi = node["metadata"]["rsmi"]
                print(f"Analyzing final step reaction: {rsmi}")

                # Check if this is any type of amide formation reaction using the checker
                amide_formation_reactions = [
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_OS",
                    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                    "Acyl chloride with ammonia to amide",
                    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                    "Acyl chloride with secondary amine to amide",
                    "Carboxylic acid with primary amine to amide",
                    "Ester with ammonia to amide",
                    "Ester with primary amine to amide",
                    "Ester with secondary amine to amide",
                    "Acylation of secondary amines with anhydrides",
                    "Acylation of secondary amines",
                    "Acylation of primary amines",
                    "Schotten-Baumann to ester",
                    "Schotten-Baumann_amide",
                    "Carboxylic acid to amide conversion",
                    "Aminolysis of esters",
                ]

                for reaction_type in amide_formation_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        print(f"Detected amide formation as final step: {reaction_type}")
                        final_step_is_amide_formation = True
                        # Don't return, continue checking other criteria

                # If no specific reaction type matched, check for reactants and products
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Checking product for amide: {product}")
                # Check if product contains amide
                product_has_amide = (
                    checker.check_fg("Primary amide", product)
                    or checker.check_fg("Secondary amide", product)
                    or checker.check_fg("Tertiary amide", product)
                )

                if product_has_amide:
                    print("Product contains amide functional group")

                    # Check if any reactant contains amide to see if it's newly formed
                    reactants_have_amide = False
                    for reactant in reactants:
                        if (
                            checker.check_fg("Primary amide", reactant)
                            or checker.check_fg("Secondary amide", reactant)
                            or checker.check_fg("Tertiary amide", reactant)
                        ):
                            reactants_have_amide = True
                            print(f"Reactant already contains amide: {reactant}")
                            break

                    # Only consider it amide formation if the amide is not present in reactants
                    if not reactants_have_amide:
                        # Check for potential amide-forming reactants
                        has_amide_precursor = False
                        for reactant in reactants:
                            if (
                                checker.check_fg("Carboxylic acid", reactant)
                                or checker.check_fg("Acyl halide", reactant)
                                or checker.check_fg("Anhydride", reactant)
                                or checker.check_fg("Ester", reactant)
                            ):
                                has_amide_precursor = True
                                print(f"Found amide precursor: {reactant}")
                                break

                        # Check for amine reactants
                        has_amine = False
                        for reactant in reactants:
                            if (
                                checker.check_fg("Primary amine", reactant)
                                or checker.check_fg("Secondary amine", reactant)
                                or checker.check_fg("Aniline", reactant)
                            ):
                                has_amine = True
                                print(f"Found amine reactant: {reactant}")
                                break

                        if has_amide_precursor and has_amine:
                            print("Detected amide formation as final step (by functional groups)")
                            final_step_is_amide_formation = True

            except Exception as e:
                print(f"Error analyzing reaction: {e}")

        # Traverse children with incremented depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    print(f"Final result: {final_step_is_amide_formation}")
    return final_step_is_amide_formation
