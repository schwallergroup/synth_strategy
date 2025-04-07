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
    Detects if the synthesis route involves amide formation in the final step.
    """
    final_amide_formation = False
    target_mol_smiles = route["smiles"]

    def dfs_traverse(node, depth=0, path=None):
        nonlocal final_amide_formation, target_mol_smiles
        if path is None:
            path = []

        # Add current node to path
        path.append(node)

        if node["type"] == "reaction" and node.get("metadata", {}).get("rsmi"):
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0]
            product_smiles = rsmi.split(">")[-1]

            # Check if this reaction produces the target molecule
            is_final_step = False
            if (
                len(path) >= 2
                and path[-2]["type"] == "mol"
                and path[-2]["smiles"] == target_mol_smiles
            ):
                is_final_step = True
                print(f"Found reaction leading to target molecule at depth {depth}")

            # Check for amide formation reactions
            is_amide_reaction = False

            # Check for specific amide formation reactions from the provided list
            amide_reaction_types = [
                "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                "Ester with primary amine to amide",
                "Ester with secondary amine to amide",
                "Ester with ammonia to amide",
                "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                "Acyl chloride with secondary amine to amide",
                "Acyl chloride with ammonia to amide",
                "Carboxylic acid with primary amine to amide",
                "Schotten-Baumann to ester",
                "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                "Acylation of primary amines",
                "Acylation of secondary amines",
            ]

            for reaction_type in amide_reaction_types:
                if checker.check_reaction(reaction_type, rsmi):
                    print(f"Detected amide formation reaction: {reaction_type}")
                    is_amide_reaction = True
                    break

            # If no specific reaction type matched, check for functional group changes
            if not is_amide_reaction:
                # Check for reactant functional groups that could form amides
                reactant_list = reactants_smiles.split(".")
                acid_present = any(
                    checker.check_fg("Carboxylic acid", r) for r in reactant_list
                )
                ester_present = any(checker.check_fg("Ester", r) for r in reactant_list)
                acyl_halide_present = any(
                    checker.check_fg("Acyl halide", r) for r in reactant_list
                )
                anhydride_present = any(
                    checker.check_fg("Anhydride", r) for r in reactant_list
                )

                # Check for amine reactants
                primary_amine_present = any(
                    checker.check_fg("Primary amine", r) for r in reactant_list
                )
                secondary_amine_present = any(
                    checker.check_fg("Secondary amine", r) for r in reactant_list
                )
                amine_present = primary_amine_present or secondary_amine_present

                # Check for amide in product
                amide_in_product = (
                    checker.check_fg("Primary amide", product_smiles)
                    or checker.check_fg("Secondary amide", product_smiles)
                    or checker.check_fg("Tertiary amide", product_smiles)
                )

                # Check for amide in reactants (to verify it's newly formed)
                amide_in_reactants = any(
                    checker.check_fg("Primary amide", r)
                    or checker.check_fg("Secondary amide", r)
                    or checker.check_fg("Tertiary amide", r)
                    for r in reactant_list
                )

                # Verify amide is formed in the reaction
                if (
                    (
                        acid_present
                        or ester_present
                        or acyl_halide_present
                        or anhydride_present
                    )
                    and amine_present
                    and amide_in_product
                    and not amide_in_reactants
                ):
                    print(f"Detected amide formation from functional group analysis")
                    is_amide_reaction = True

            # If this is an amide formation reaction and it's the final step
            if is_amide_reaction and is_final_step:
                print(f"Detected late-stage amide formation at depth {depth}")
                final_amide_formation = True

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1, path.copy())

    dfs_traverse(route)
    print(f"Final result: late-stage amide formation = {final_amide_formation}")
    return final_amide_formation
