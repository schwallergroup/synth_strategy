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
    This function detects if the final step in the synthesis is an amide formation.
    """
    strategy_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal strategy_detected

        # Check if this is a late-stage reaction (depth 0, 1, or 2)
        if node["type"] == "reaction" and depth <= 2:
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if product contains amide
                if (
                    checker.check_fg("Primary amide", product_smiles)
                    or checker.check_fg("Secondary amide", product_smiles)
                    or checker.check_fg("Tertiary amide", product_smiles)
                ):

                    # Check if this is an amide formation reaction
                    amide_formation_reactions = [
                        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                        "Acyl chloride with ammonia to amide",
                        "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                        "Acyl chloride with secondary amine to amide",
                        "Carboxylic acid with primary amine to amide",
                        "Ester with ammonia to amide",
                        "Ester with primary amine to amide",
                        "Ester with secondary amine to amide",
                        "Schotten-Baumann_amide",
                        "Nitrile to amide",
                    ]

                    for reaction_type in amide_formation_reactions:
                        if checker.check_reaction(reaction_type, rsmi):
                            print(
                                f"Detected late-stage amide formation via {reaction_type} at depth {depth}"
                            )
                            strategy_detected = True
                            return

                    # If no specific reaction type matched, check for reactants with acid/amine
                    has_acid = False
                    has_amine = False
                    has_acyl_halide = False
                    has_ester = False

                    for reactant_smiles in reactants_smiles:
                        if checker.check_fg("Carboxylic acid", reactant_smiles):
                            has_acid = True
                        if checker.check_fg(
                            "Primary amine", reactant_smiles
                        ) or checker.check_fg("Secondary amine", reactant_smiles):
                            has_amine = True
                        if checker.check_fg("Acyl halide", reactant_smiles):
                            has_acyl_halide = True
                        if checker.check_fg("Ester", reactant_smiles):
                            has_ester = True

                    # Check for common amide formation patterns
                    if (
                        (has_acid and has_amine)
                        or (has_acyl_halide and has_amine)
                        or (has_ester and has_amine)
                    ):
                        # Verify amide is newly formed by checking reactants don't already have amide
                        reactants_have_amide = False
                        for reactant_smiles in reactants_smiles:
                            if (
                                checker.check_fg("Primary amide", reactant_smiles)
                                or checker.check_fg("Secondary amide", reactant_smiles)
                                or checker.check_fg("Tertiary amide", reactant_smiles)
                            ):
                                reactants_have_amide = True
                                break

                        if not reactants_have_amide:
                            print(
                                f"Detected late-stage amide formation from functional groups at depth {depth}"
                            )
                            strategy_detected = True
                            return

        # Traverse children with increased depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)
    return strategy_detected
