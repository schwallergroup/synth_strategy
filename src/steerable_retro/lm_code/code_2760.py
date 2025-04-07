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
    This function detects if the synthetic route involves protection of a carboxylic acid as an ester
    and subsequent deprotection.
    """
    # Track both protection and deprotection events
    protection_events = []
    deprotection_events = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check for protection reactions (acid → ester in forward direction)
                if checker.check_reaction("Protection of carboxylic acid", rsmi):
                    print(f"Found carboxylic acid protection reaction at depth {depth}: {rsmi}")
                    protection_events.append((depth, rsmi))
                    return

                # Check for deprotection reactions (ester → acid in forward direction)
                if (
                    checker.check_reaction("Ester saponification (methyl deprotection)", rsmi)
                    or checker.check_reaction("Ester saponification (alkyl deprotection)", rsmi)
                    or checker.check_reaction("COOH ethyl deprotection", rsmi)
                    or checker.check_reaction("Deprotection of carboxylic acid", rsmi)
                ):
                    print(f"Found carboxylic acid deprotection reaction at depth {depth}: {rsmi}")
                    deprotection_events.append((depth, rsmi))
                    return

                # Alternative check for protection: carboxylic acid to ester transformation
                for reactant in reactants_smiles:
                    if checker.check_fg("Carboxylic acid", reactant):
                        product = product_smiles
                        if checker.check_fg("Ester", product):
                            print(
                                f"Found carboxylic acid to ester transformation at depth {depth}: {rsmi}"
                            )
                            protection_events.append((depth, rsmi))
                            return

                # Alternative check for deprotection: ester to carboxylic acid transformation
                for reactant in reactants_smiles:
                    if checker.check_fg("Ester", reactant) and checker.check_fg(
                        "Carboxylic acid", product_smiles
                    ):
                        print(
                            f"Found ester to carboxylic acid transformation at depth {depth}: {rsmi}"
                        )
                        deprotection_events.append((depth, rsmi))
                        return

            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # A protection strategy requires both protection and deprotection events
    has_protection_strategy = len(protection_events) > 0 or len(deprotection_events) > 0

    print(f"Protection events: {len(protection_events)}")
    print(f"Deprotection events: {len(deprotection_events)}")
    print(f"Carboxylic acid protection strategy detected: {has_protection_strategy}")

    return has_protection_strategy
