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
    Detects if the route contains a reductive amination step (primary amine + aldehyde -> secondary amine).
    """
    has_reductive_amination = False

    def dfs_traverse(node):
        nonlocal has_reductive_amination

        if node["type"] == "reaction" and not has_reductive_amination:
            try:
                rsmi = node["metadata"]["rsmi"]

                # Direct check for reductive amination reaction
                if checker.check_reaction("Reductive amination with aldehyde", rsmi):
                    print(f"Found reductive amination via reaction check: {rsmi}")
                    has_reductive_amination = True
                    return

                # If direct check fails, perform more detailed analysis
                reactants_part = rsmi.split(">")[0]
                product_part = rsmi.split(">")[-1]
                reactant_fragments = reactants_part.split(".")

                # Check if we have multiple reactants
                if len(reactant_fragments) > 1:
                    # Check for primary amine in reactants
                    primary_amine_fragments = [
                        frag
                        for frag in reactant_fragments
                        if checker.check_fg("Primary amine", frag)
                    ]

                    # Check for aldehyde in reactants
                    aldehyde_fragments = [
                        frag for frag in reactant_fragments if checker.check_fg("Aldehyde", frag)
                    ]

                    # Check for secondary amine in product
                    has_secondary_amine = checker.check_fg("Secondary amine", product_part)

                    # Verify we have both primary amine and aldehyde in different fragments
                    # and the product contains a secondary amine
                    if (
                        primary_amine_fragments
                        and aldehyde_fragments
                        and has_secondary_amine
                        and set(primary_amine_fragments) != set(aldehyde_fragments)
                    ):

                        print(f"Found reductive amination via fragment analysis: {rsmi}")
                        print(f"Primary amine fragments: {primary_amine_fragments}")
                        print(f"Aldehyde fragments: {aldehyde_fragments}")
                        print(f"Product with secondary amine: {product_part}")
                        has_reductive_amination = True
            except Exception as e:
                print(f"Error analyzing reaction: {e}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    print(f"Reductive amination strategy detected: {has_reductive_amination}")
    return has_reductive_amination
