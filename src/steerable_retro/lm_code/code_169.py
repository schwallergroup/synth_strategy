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
    Detects if the synthetic route involves multiple esterification reactions.
    """
    esterification_count = 0
    processed_reactions = set()  # To avoid counting the same reaction twice

    def dfs_traverse(node, depth=0):
        nonlocal esterification_count

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]

            # Skip if we've already processed this reaction
            if rsmi in processed_reactions:
                return
            processed_reactions.add(rsmi)

            # Check for various esterification reaction types
            esterification_reactions = [
                "Esterification of Carboxylic Acids",
                "Schotten-Baumann to ester",
                "O-alkylation of carboxylic acids with diazo compounds",
                "Transesterification",
                "Oxidative esterification of primary alcohols",
                "Acetic anhydride and alcohol to ester",
                "Mitsunobu esterification",
                "Pinner reaction to ester",
                "Oxidation of alcohol and aldehyde to ester",
            ]

            # Check for esterification formation reactions
            for reaction_type in esterification_reactions:
                if checker.check_reaction(reaction_type, rsmi):
                    print(
                        f"Detected esterification reaction ({reaction_type}) at depth {depth}: {rsmi}"
                    )
                    esterification_count += 1
                    break  # Count each reaction only once

            # Check for ester hydrolysis reactions
            hydrolysis_reactions = [
                "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
                "Ester saponification (methyl deprotection)",
                "Ester saponification (alkyl deprotection)",
                "COOH ethyl deprotection",
            ]

            for reaction_type in hydrolysis_reactions:
                if checker.check_reaction(reaction_type, rsmi):
                    print(
                        f"Detected ester hydrolysis reaction ({reaction_type}) at depth {depth}: {rsmi}"
                    )
                    esterification_count += 1
                    break  # Count each reaction only once

            # If no specific reaction type matched, use fallback method
            if all(
                not checker.check_reaction(r, rsmi)
                for r in esterification_reactions + hydrolysis_reactions
            ):
                reactants_str = rsmi.split(">")[0]
                product_str = rsmi.split(">")[-1]
                reactants = reactants_str.split(".")

                # Check for ester formation
                has_ester_product = checker.check_fg("Ester", product_str)
                has_ester_reactant = any(checker.check_fg("Ester", r) for r in reactants)

                # Check for ester formation
                if has_ester_product and not has_ester_reactant:
                    # Check various esterification pathways

                    # 1. Carboxylic acid + alcohol → ester
                    has_carboxylic_acid = any(
                        checker.check_fg("Carboxylic acid", r) for r in reactants
                    )
                    has_alcohol = any(
                        checker.check_fg("Primary alcohol", r)
                        or checker.check_fg("Secondary alcohol", r)
                        or checker.check_fg("Tertiary alcohol", r)
                        or checker.check_fg("Aromatic alcohol", r)
                        or checker.check_fg("Phenol", r)
                        for r in reactants
                    )

                    # 2. Acyl halide + alcohol → ester
                    has_acyl_halide = any(checker.check_fg("Acyl halide", r) for r in reactants)

                    # 3. Anhydride + alcohol → ester
                    has_anhydride = any(checker.check_fg("Anhydride", r) for r in reactants)

                    # 4. Aldehyde oxidation to ester
                    has_aldehyde = any(checker.check_fg("Aldehyde", r) for r in reactants)

                    if (
                        (has_carboxylic_acid and has_alcohol)
                        or (has_acyl_halide and has_alcohol)
                        or (has_anhydride and has_alcohol)
                        or has_aldehyde
                    ):
                        print(
                            f"Detected esterification reaction (fallback method) at depth {depth}: {rsmi}"
                        )
                        esterification_count += 1

                # Check for ester hydrolysis/consumption
                elif has_ester_reactant and not has_ester_product:
                    has_carboxylic_acid_product = checker.check_fg("Carboxylic acid", product_str)
                    has_alcohol_product = any(
                        checker.check_fg("Primary alcohol", product_str)
                        or checker.check_fg("Secondary alcohol", product_str)
                        or checker.check_fg("Tertiary alcohol", product_str)
                        or checker.check_fg("Aromatic alcohol", product_str)
                        or checker.check_fg("Phenol", product_str)
                        for r in reactants
                    )

                    if has_carboxylic_acid_product or has_alcohol_product:
                        print(
                            f"Detected ester hydrolysis (fallback method) at depth {depth}: {rsmi}"
                        )
                        esterification_count += 1

                # Check for transesterification
                elif has_ester_reactant and has_ester_product:
                    has_alcohol = any(
                        checker.check_fg("Primary alcohol", r)
                        or checker.check_fg("Secondary alcohol", r)
                        or checker.check_fg("Tertiary alcohol", r)
                        or checker.check_fg("Aromatic alcohol", r)
                        or checker.check_fg("Phenol", r)
                        for r in reactants
                    )

                    if has_alcohol:
                        print(
                            f"Detected transesterification (fallback method) at depth {depth}: {rsmi}"
                        )
                        esterification_count += 1

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    print(f"Total esterification reactions found: {esterification_count}")
    return esterification_count >= 2
