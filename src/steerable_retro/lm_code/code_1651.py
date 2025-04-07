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

root_data = "/home/andres/Documents/steerable_retro/data"

fg_args = {
    "file_path": f"{root_data}/patterns/functional_groups.json",
    "value_field": "pattern",
    "key_field": "name",
}
reaction_class_args = {
    "file_path": f"{root_data}/patterns/smirks.json",
    "value_field": "smirks",
    "key_field": "name",
}
ring_smiles_args = {
    "file_path": f"{root_data}/patterns/chemical_rings_smiles.json",
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
    This function detects if the synthesis route employs a late-stage amination strategy,
    where an NH2 group is introduced in the final step.
    """
    late_stage_amination_found = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_amination_found

        # Check if this is a reaction node
        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            # Late stage means depth 0, 1, or 2 (close to target molecule)
            if depth <= 2:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Examining reaction at depth {depth}: {rsmi}")

                # Check if primary amine is in the product
                if checker.check_fg("Primary amine", product):
                    print(f"Found primary amine in product at depth {depth}")

                    # Check if primary amine is NOT in any of the reactants
                    primary_amine_in_reactants = any(
                        checker.check_fg("Primary amine", reactant) for reactant in reactants
                    )

                    if not primary_amine_in_reactants:
                        print(f"Primary amine not found in reactants at depth {depth}")

                        # Verify this is an amination reaction
                        amination_reactions = [
                            "Reduction of nitro groups to amines",
                            "Reductive amination with aldehyde",
                            "Reductive amination with ketone",
                            "Reductive amination with alcohol",
                            "N-alkylation of primary amines with alkyl halides",
                            "Reduction of nitrile to amine",
                            "Reduction of primary amides to amines",
                            "Azide to amine reduction (Staudinger)",
                            "Phthalimide deprotection",
                            "N-glutarimide deprotection",
                            "Boc amine deprotection",
                            "Boc amine deprotection to NH-NH2",
                            "Tert-butyl deprotection of amine",
                            "Hydrazine synthesis from amine",
                            "Amine to azide",  # This can be the reverse direction
                            "Primary amine to fluoride",  # This can be the reverse direction
                            "Primary amine to chloride",  # This can be the reverse direction
                            "Primary amine to bromide",  # This can be the reverse direction
                            "Primary amine to iodide",  # This can be the reverse direction
                            "Aminolysis of esters",
                            "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                            "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                            "Acylation of primary amines",
                            "Acylation of secondary amines",
                            "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
                            "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
                            "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
                            "Goldberg coupling aryl amine-aryl chloride",
                            "Ullmann-Goldberg Substitution amine",
                            "Ring opening of epoxide with amine",
                            "Hydrogenolysis of amides/imides/carbamates",
                            "Hydrolysis of amides/imides/carbamates",
                        ]

                        for reaction_type in amination_reactions:
                            if checker.check_reaction(reaction_type, rsmi):
                                print(
                                    f"Late-stage amination detected: {reaction_type} at depth {depth}"
                                )
                                late_stage_amination_found = True
                                break

                        # Additional check for nitro reduction
                        if not late_stage_amination_found:
                            if any(
                                checker.check_fg("Nitro group", reactant) for reactant in reactants
                            ) and not checker.check_fg("Nitro group", product):
                                print(
                                    f"Late-stage nitro reduction to amine detected at depth {depth}"
                                )
                                late_stage_amination_found = True

                        # Additional check for nitrile reduction
                        if not late_stage_amination_found:
                            if any(
                                checker.check_fg("Nitrile", reactant) for reactant in reactants
                            ) and not checker.check_fg("Nitrile", product):
                                print(
                                    f"Late-stage nitrile reduction to amine detected at depth {depth}"
                                )
                                late_stage_amination_found = True

                        # Additional check for azide reduction
                        if not late_stage_amination_found:
                            if any(
                                checker.check_fg("Azide", reactant) for reactant in reactants
                            ) and not checker.check_fg("Azide", product):
                                print(
                                    f"Late-stage azide reduction to amine detected at depth {depth}"
                                )
                                late_stage_amination_found = True

                        # Fallback: If we've identified a primary amine in the product but not in reactants,
                        # and we're at a late stage (depth ≤ 2), consider this a late-stage amination
                        if not late_stage_amination_found:
                            print(f"Unclassified late-stage amination detected at depth {depth}")
                            late_stage_amination_found = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    return late_stage_amination_found
