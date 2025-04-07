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
    This function detects a synthesis strategy where the final step is an amide coupling.
    """
    final_step_is_amide_coupling = False

    def dfs_traverse(node, depth=0):
        nonlocal final_step_is_amide_coupling

        # Check if this is a reaction node at depth 0 (final reaction)
        if node["type"] == "reaction" and depth == 0:
            print(f"Examining final reaction at depth {depth}")
            # Check if this is an amide coupling reaction
            rsmi = node["metadata"]["rsmi"]
            print(f"Reaction SMILES: {rsmi}")

            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]
            print(f"Product: {product}")
            print(f"Reactants: {reactants}")

            # Check if this is a known amide coupling reaction type
            amide_coupling_reactions = [
                "Carboxylic acid with primary amine to amide",
                "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                "Ester with primary amine to amide",
                "Ester with secondary amine to amide",
                "Acyl chloride with secondary amine to amide",
                "Schotten-Baumann_amide",
            ]

            for reaction_type in amide_coupling_reactions:
                if checker.check_reaction(reaction_type, rsmi):
                    print(f"Detected amide coupling reaction: {reaction_type}")
                    final_step_is_amide_coupling = True
                    return

            # If no specific reaction type matched, check for functional groups
            # Check for amide in product
            amide_types = ["Primary amide", "Secondary amide", "Tertiary amide"]
            has_amide = any(checker.check_fg(amide_type, product) for amide_type in amide_types)

            if has_amide:
                print("Product contains amide group")

                # Check for carboxylic acid or acyl halide in reactants
                acid_types = ["Carboxylic acid", "Acyl halide"]
                has_acid_or_acyl = any(
                    any(checker.check_fg(acid_type, r) for acid_type in acid_types)
                    for r in reactants
                )

                # Check for amine in reactants
                amine_types = ["Primary amine", "Secondary amine"]
                has_amine = any(
                    any(checker.check_fg(amine_type, r) for amine_type in amine_types)
                    for r in reactants
                )

                print(f"Reactants contain acid/acyl: {has_acid_or_acyl}, amine: {has_amine}")

                if has_acid_or_acyl and has_amine:
                    print("Detected amide coupling based on functional groups")
                    final_step_is_amide_coupling = True

        # If this is the root molecule node, check its child reactions
        if node["type"] == "mol" and depth == 0 and "children" in node:
            print("Found root molecule node, checking its child reactions")
            for child in node.get("children", []):
                # Keep depth at 0 for immediate child reactions of the root
                dfs_traverse(child, 0)
        else:
            # Continue normal traversal for other nodes
            for child in node.get("children", []):
                dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    print(f"Final result: {final_step_is_amide_coupling}")
    return final_step_is_amide_coupling
