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
    This function detects if the synthesis uses a late-stage amide coupling strategy,
    specifically looking for amide formation in the final step.
    """
    amide_formation_at_final_step = False

    def dfs_traverse(node, current_depth=0):
        nonlocal amide_formation_at_final_step

        print(f"Traversing node at depth {current_depth}, type: {node['type']}")

        if node["type"] == "reaction" and current_depth == 1:
            print("Analyzing final step reaction (depth 1)")
            # Check if this is a depth 1 reaction (final step)
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Final reaction SMILES: {rsmi}")
                print(f"Reactants: {reactants}")
                print(f"Product: {product}")

                # Check for amide coupling reactions
                amide_coupling_reactions = [
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_OS",
                    "Acyl chloride with ammonia to amide",
                    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                    "Acyl chloride with primary amine to imide",
                    "Acyl chloride with secondary amine to amide",
                    "Carboxylic acid with primary amine to amide",
                    "Ester with ammonia to amide",
                    "Ester with primary amine to amide",
                    "Ester with secondary amine to amide",
                    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                    "Schotten-Baumann_amide",
                    "Acylation of primary amines",
                    "Acylation of secondary amines",
                    "Carboxylic acid to amide conversion",
                    "Nitrile and hydrogen peroxide to amide",
                ]

                # Check if any of the amide coupling reactions match
                for reaction_type in amide_coupling_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        print(f"Detected {reaction_type} in final step")
                        amide_formation_at_final_step = True
                        return

                # If no specific reaction matched, check for amide formation by comparing functional groups
                try:
                    # Check if product has amide
                    amide_in_product = (
                        checker.check_fg("Primary amide", product)
                        or checker.check_fg("Secondary amide", product)
                        or checker.check_fg("Tertiary amide", product)
                    )

                    print(f"Product contains amide: {amide_in_product}")

                    if amide_in_product:
                        # Check if reactants have necessary precursors for amide formation
                        has_amine = False
                        has_carboxylic_acid_or_derivative = False

                        for reactant in reactants:
                            # Check for amine precursors
                            if (
                                checker.check_fg("Primary amine", reactant)
                                or checker.check_fg("Secondary amine", reactant)
                                or checker.check_fg("Aniline", reactant)
                            ):
                                has_amine = True
                                print(f"Found amine in reactant: {reactant}")

                            # Check for carboxylic acid or derivatives
                            if (
                                checker.check_fg("Carboxylic acid", reactant)
                                or checker.check_fg("Acyl halide", reactant)
                                or checker.check_fg("Ester", reactant)
                                or checker.check_fg("Anhydride", reactant)
                            ):
                                has_carboxylic_acid_or_derivative = True
                                print(
                                    f"Found carboxylic acid or derivative in reactant: {reactant}"
                                )

                        # Check if at least one reactant doesn't have the amide (new amide formation)
                        new_amide_formation = False
                        for reactant in reactants:
                            reactant_has_amide = (
                                checker.check_fg("Primary amide", reactant)
                                or checker.check_fg("Secondary amide", reactant)
                                or checker.check_fg("Tertiary amide", reactant)
                            )
                            if not reactant_has_amide:
                                new_amide_formation = True
                                break

                        if new_amide_formation and has_amine and has_carboxylic_acid_or_derivative:
                            print("Detected amide formation in final step (FG analysis)")
                            amide_formation_at_final_step = True
                            return
                except Exception as e:
                    print(f"Error in amide functional group detection: {e}")

        # Process children with increased depth
        for child in node.get("children", []):
            dfs_traverse(child, current_depth + 1)

    # Start traversal
    dfs_traverse(route)
    print(f"Final result: {amide_formation_at_final_step}")
    return amide_formation_at_final_step
