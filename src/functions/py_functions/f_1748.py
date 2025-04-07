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
    Detects reduction of an ester to a primary alcohol.
    """
    found_ester_reduction = False

    def dfs_traverse(node):
        nonlocal found_ester_reduction

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            products_part = rsmi.split(">")[-1]

            print(f"Examining reaction: {rsmi}")

            # First, check if this is a known ester reduction reaction
            if checker.check_reaction("Reduction of ester to primary alcohol", rsmi):
                print("Found ester reduction reaction via reaction check")
                found_ester_reduction = True
                return

            # Check for functional group transformation
            if checker.check_fg("Ester", reactants_part) and checker.check_fg(
                "Primary alcohol", products_part
            ):

                # Make sure we're not just seeing an ester and alcohol in different parts
                if not checker.check_fg(
                    "Primary alcohol", reactants_part
                ) or not checker.check_fg("Ester", products_part):
                    print(
                        "Found potential ester reduction: ester in reactants and primary alcohol in products"
                    )

                    # Check for reduction reactions that could convert esters to alcohols
                    reduction_reactions = [
                        "Reduction of ester to primary alcohol",
                        "Reduction of carboxylic acid to primary alcohol",
                        "Reduction of aldehydes and ketones to alcohols",  # This might catch some related reductions
                    ]

                    for rxn_type in reduction_reactions:
                        if checker.check_reaction(rxn_type, rsmi):
                            print(f"Confirmed ester reduction via {rxn_type}")
                            found_ester_reduction = True
                            return

                    # If no specific reaction type is found, look at the reaction more carefully
                    # Check if the reaction involves a reducing agent
                    reducing_agents = [
                        "[Al]",
                        "[H]",
                        "BH4",
                        "AlH4",
                        "LiAlH4",
                        "NaBH4",
                        "DIBAL",
                        "LAH",
                    ]
                    reagents_part = rsmi.split(">")[1]

                    for agent in reducing_agents:
                        if agent in reagents_part:
                            print(f"Found reducing agent {agent} in reaction")
                            found_ester_reduction = True
                            return

                    # As a last resort, check if this is a general reduction
                    # by looking for the pattern of ester being converted to alcohol
                    # without other major functional group changes
                    if not checker.check_fg(
                        "Carboxylic acid", reactants_part
                    ) and not checker.check_fg("Aldehyde", reactants_part):
                        print(
                            "Detected ester to alcohol transformation without other major functional group changes"
                        )
                        found_ester_reduction = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return found_ester_reduction
