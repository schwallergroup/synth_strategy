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
    Detects if the synthesis route uses alcohol activation (e.g., mesylation, tosylation, triflation, halogenation)
    followed by nucleophilic substitution.
    """
    activation_reactions = []
    substitution_reactions = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            print(f"Examining reaction at depth {depth}: {rsmi}")
            print(f"Reactants: {reactants}")
            print(f"Product: {product}")

            # Check for alcohol activation reactions
            has_alcohol = any(
                checker.check_fg("Primary alcohol", r)
                or checker.check_fg("Secondary alcohol", r)
                or checker.check_fg("Tertiary alcohol", r)
                or checker.check_fg("Aromatic alcohol", r)
                for r in reactants
            )

            if has_alcohol:
                print(f"Found alcohol in reactants at depth {depth}")

            has_activated = (
                checker.check_fg("Mesylate", product)
                or checker.check_fg("Tosylate", product)
                or checker.check_fg("Triflate", product)
                or checker.check_fg("Primary halide", product)
                or checker.check_fg("Secondary halide", product)
            )

            if has_activated:
                print(f"Found activated group in product at depth {depth}")

            # Check for specific activation reactions if available
            is_activation_rxn = (
                checker.check_reaction("Formation of Sulfonic Esters", rsmi)
                or checker.check_reaction(
                    "Formation of Sulfonic Esters on TMS protected alcohol", rsmi
                )
                or checker.check_reaction("Alcohol to triflate conversion", rsmi)
                or checker.check_reaction("Alcohol to chloride_sulfonyl chloride", rsmi)
                or checker.check_reaction("Alcohol to chloride_SOCl2", rsmi)
                or checker.check_reaction("Alcohol to chloride_PCl5_ortho", rsmi)
                or checker.check_reaction("Alcohol to chloride_POCl3", rsmi)
                or checker.check_reaction("Alcohol to chloride_HCl", rsmi)
                or checker.check_reaction("Appel reaction", rsmi)
            )

            if is_activation_rxn:
                print(f"Found activation reaction at depth {depth}")

            # Detect activation if either pattern is found
            if (has_alcohol and has_activated) or is_activation_rxn:
                activation_reactions.append((depth, rsmi))
                print(f"Found alcohol activation at depth {depth}: {rsmi}")

            # Check for nucleophilic substitution with activated leaving groups
            has_leaving_group = any(
                checker.check_fg("Mesylate", r)
                or checker.check_fg("Tosylate", r)
                or checker.check_fg("Triflate", r)
                or checker.check_fg("Primary halide", r)
                or checker.check_fg("Secondary halide", r)
                for r in reactants
            )

            if has_leaving_group:
                print(f"Found leaving group in reactants at depth {depth}")

            # Check if the leaving group is gone in the product
            leaving_group_gone = not (
                checker.check_fg("Mesylate", product)
                or checker.check_fg("Tosylate", product)
                or checker.check_fg("Triflate", product)
                or checker.check_fg("Primary halide", product)
                or checker.check_fg("Secondary halide", product)
            )

            if leaving_group_gone:
                print(f"Leaving group is gone in product at depth {depth}")

            # Check for nucleophilic substitution reactions
            is_subst_rxn = (
                checker.check_reaction("Williamson Ether Synthesis", rsmi)
                or checker.check_reaction("S-alkylation of thiols", rsmi)
                or checker.check_reaction("N-alkylation of primary amines with alkyl halides", rsmi)
                or checker.check_reaction(
                    "N-alkylation of secondary amines with alkyl halides", rsmi
                )
                or checker.check_reaction("Mitsunobu aryl ether", rsmi)
                or checker.check_reaction("Mitsunobu esterification", rsmi)
                or checker.check_reaction("Alcohol to azide", rsmi)
            )

            if is_subst_rxn:
                print(f"Found substitution reaction at depth {depth}")

            # Detect substitution if either pattern is found
            if (has_leaving_group and leaving_group_gone) or is_subst_rxn:
                substitution_reactions.append((depth, rsmi))
                print(f"Found nucleophilic substitution at depth {depth}: {rsmi}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    print(f"Activation reactions found: {len(activation_reactions)}")
    print(f"Substitution reactions found: {len(substitution_reactions)}")

    # Check if we found both activation and substitution reactions
    if activation_reactions and substitution_reactions:
        # Sort by depth
        activation_reactions.sort(key=lambda x: x[0])
        substitution_reactions.sort(key=lambda x: x[0])

        print(f"Sorted activation reactions: {activation_reactions}")
        print(f"Sorted substitution reactions: {substitution_reactions}")

        # In forward synthesis, activation should happen before substitution
        # In retrosynthesis (our traversal), activation depth should be greater than substitution depth
        for act_depth, act_rsmi in activation_reactions:
            for subst_depth, subst_rsmi in substitution_reactions:
                if act_depth > subst_depth:  # Changed from < to >
                    print(f"Alcohol activation-substitution strategy detected:")
                    print(f"  Activation at depth {act_depth}: {act_rsmi}")
                    print(f"  Substitution at depth {subst_depth}: {subst_rsmi}")
                    return True

    return False
