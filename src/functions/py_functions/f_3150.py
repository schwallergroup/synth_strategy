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
    Detects the use of ester formation as a protection strategy for carboxylic acids.
    """
    # Track protection and deprotection events
    protection_events = []
    deprotection_events = []

    def has_other_reactive_groups(smiles):
        """Check if molecule has other reactive functional groups that might need protection"""
        reactive_groups = [
            "Primary alcohol",
            "Secondary alcohol",
            "Tertiary alcohol",
            "Primary amine",
            "Secondary amine",
            "Aldehyde",
            "Ketone",
            "Thiol",
            "Phenol",
            "Alkyne",
            "Alkene",
        ]
        for group in reactive_groups:
            if checker.check_fg(group, smiles):
                print(f"Found reactive group: {group} in {smiles}")
                return True
        return False

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for esterification (protection)
                is_esterification = checker.check_reaction(
                    "Esterification of Carboxylic Acids", rsmi
                )

                if is_esterification:
                    # Verify carboxylic acid in reactants and ester in product
                    for reactant in reactants:
                        if checker.check_fg(
                            "Carboxylic acid", reactant
                        ) and checker.check_fg("Ester", product):
                            print(
                                f"Detected ester formation (potential protection) at depth {depth}"
                            )
                            # Store the depth and the molecules involved
                            protection_events.append((depth, reactant, product))

                # Check for deprotection
                is_deprotection = (
                    checker.check_reaction(
                        "Ester saponification (methyl deprotection)", rsmi
                    )
                    or checker.check_reaction(
                        "Ester saponification (alkyl deprotection)", rsmi
                    )
                    or checker.check_reaction(
                        "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
                        rsmi,
                    )
                )

                if is_deprotection:
                    # Verify ester in reactants and carboxylic acid in product
                    for reactant in reactants:
                        if checker.check_fg("Ester", reactant) and checker.check_fg(
                            "Carboxylic acid", product
                        ):
                            print(f"Detected ester deprotection at depth {depth}")
                            # Store the depth and the molecules involved
                            deprotection_events.append((depth, reactant, product))

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Analyze if protection strategy was used
    print(
        f"Found {len(protection_events)} protection events and {len(deprotection_events)} deprotection events"
    )

    # Case 1: Both protection and deprotection events exist
    if protection_events and deprotection_events:
        # In retrosynthesis, protection should have higher depth than deprotection
        # (protection happens earlier in forward synthesis)
        for prot_depth, prot_reactant, prot_product in protection_events:
            for deprot_depth, deprot_reactant, deprot_product in deprotection_events:
                if prot_depth > deprot_depth:
                    print(
                        f"Confirmed ester formation used as protection strategy: protection at depth {prot_depth}, deprotection at depth {deprot_depth}"
                    )
                    return True

    # Case 2: Only protection events exist
    if protection_events:
        # Check if these protection events are likely for protection rather than the main goal
        for depth, reactant, product in protection_events:
            if has_other_reactive_groups(reactant) or has_other_reactive_groups(
                product
            ):
                print(
                    f"Detected ester formation with other reactive groups present, likely for protection"
                )
                return True

        # If we have multiple protection events, it's more likely to be a protection strategy
        if len(protection_events) > 1:
            print(f"Multiple ester formations detected, likely a protection strategy")
            return True

    # Case 3: Check if there are intermediates with both ester and other reactive groups
    # This would suggest the ester is protecting a carboxylic acid while other chemistry happens
    def check_intermediates(node, depth=0, has_ester=False):
        if node["type"] == "mol" and "smiles" in node:
            mol_smiles = node["smiles"]
            if checker.check_fg("Ester", mol_smiles) and has_other_reactive_groups(
                mol_smiles
            ):
                print(
                    f"Found intermediate with ester and other reactive groups at depth {depth}"
                )
                return True

        for child in node.get("children", []):
            if check_intermediates(child, depth + 1, has_ester):
                return True

        return False

    if check_intermediates(route):
        return True

    print("No ester protection strategy detected")
    return False
