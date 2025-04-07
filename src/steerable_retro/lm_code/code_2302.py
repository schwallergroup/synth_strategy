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
    Detects if the route employs a Grignard or organolithium addition strategy
    where an aryl halide reacts with an aldehyde to form a secondary alcohol.
    """
    found_addition = False

    def dfs_traverse(node, depth=0):
        nonlocal found_addition

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]

            # Check if this is a Grignard reaction with aldehyde
            if checker.check_reaction("Grignard from aldehyde to alcohol", rsmi):
                print(f"Found Grignard addition to aldehyde at depth {depth}")
                found_addition = True
                return

            # Check for organolithium preparation followed by aldehyde addition
            if checker.check_reaction("Preparation of organolithium compounds", rsmi):
                print(f"Found organolithium preparation at depth {depth}")

                # Look for subsequent aldehyde addition in parent reaction
                for sibling in node.get("children", []):
                    if (
                        sibling["type"] == "reaction"
                        and "metadata" in sibling
                        and "rsmi" in sibling["metadata"]
                    ):
                        sibling_rsmi = sibling["metadata"]["rsmi"]
                        reactants = sibling_rsmi.split(">")[0].split(".")
                        product = sibling_rsmi.split(">")[-1]

                        # Check if any reactant is an aldehyde
                        has_aldehyde = any(checker.check_fg("Aldehyde", r) for r in reactants)

                        # Check if product is a secondary alcohol
                        if has_aldehyde and checker.check_fg("Secondary alcohol", product):
                            print(f"Found organolithium addition to aldehyde at depth {depth+1}")
                            found_addition = True
                            return

            # Manual check for aryl/alkyl halide + aldehyde → secondary alcohol
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for halide and aldehyde in reactants
            has_halide = any(
                checker.check_fg("Primary halide", r)
                or checker.check_fg("Secondary halide", r)
                or checker.check_fg("Tertiary halide", r)
                or checker.check_fg("Aromatic halide", r)
                for r in reactants
            )

            has_aldehyde = any(checker.check_fg("Aldehyde", r) for r in reactants)

            # Check for secondary alcohol in product
            has_sec_alcohol = checker.check_fg("Secondary alcohol", product)

            if has_halide and has_aldehyde and has_sec_alcohol:
                print(f"Found potential Grignard/organolithium addition at depth {depth}")
                found_addition = True

        # Continue DFS traversal
        for child in node.get("children", []):
            if not found_addition:  # Stop traversal if we already found what we're looking for
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return found_addition
