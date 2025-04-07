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
    This function detects if the synthesis uses a protection-deprotection strategy,
    looking for various protection and deprotection reactions.
    """
    # Initialize tracking variables
    protection_events = []
    deprotection_events = []

    # Define protection and deprotection reaction types to check
    protection_reactions = [
        "Boc amine protection",
        "Boc amine protection explicit",
        "Boc amine protection with Boc anhydride",
        "Boc amine protection (ethyl Boc)",
        "Boc amine protection of secondary amine",
        "Boc amine protection of primary amine",
        "Alcohol protection with silyl ethers",
    ]

    deprotection_reactions = [
        "Phthalimide deprotection",
        "Boc amine deprotection",
        "Boc amine deprotection of guanidine",
        "Boc amine deprotection to NH-NH2",
        "Tert-butyl deprotection of amine",
        "Alcohol deprotection from silyl ethers",
        "Alcohol deprotection from silyl ethers (double)",
        "Alcohol deprotection from silyl ethers (diol)",
        "N-glutarimide deprotection",
        "Hydroxyl benzyl deprotection",
        "Carboxyl benzyl deprotection",
    ]

    # Define protection and deprotection functional group pairs
    protection_fg_pairs = [
        ("Primary amine", "Phthalimide"),
        ("Primary amine", "Secondary amide"),
        ("Primary amine", "Tertiary amide"),
        ("Secondary amine", "Tertiary amide"),
        ("Primary alcohol", "Ether"),
        ("Secondary alcohol", "Ether"),
        ("Tertiary alcohol", "Ether"),
        ("Carboxylic acid", "Ester"),
    ]

    def dfs_traverse(node, depth=0):
        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for protection reactions
            for reaction in protection_reactions:
                if checker.check_reaction(reaction, rsmi):
                    print(f"Detected {reaction} at depth {depth}")
                    protection_events.append((reaction, depth))

            # Check for deprotection reactions
            for reaction in deprotection_reactions:
                if checker.check_reaction(reaction, rsmi):
                    print(f"Detected {reaction} at depth {depth}")
                    deprotection_events.append((reaction, depth))

            # If no specific reaction detected, check for functional group changes
            if not any(
                checker.check_reaction(r, rsmi)
                for r in protection_reactions + deprotection_reactions
            ):
                try:
                    product_mol = product
                    reactant_mols = reactants

                    # Check for protection by functional group changes
                    for fg_start, fg_end in protection_fg_pairs:
                        if any(
                            checker.check_fg(fg_start, r) for r in reactant_mols
                        ) and checker.check_fg(fg_end, product_mol):
                            print(
                                f"Detected potential protection: {fg_start} -> {fg_end} at depth {depth}"
                            )
                            protection_events.append((f"{fg_start} -> {fg_end}", depth))

                    # Check for deprotection by functional group changes (reverse of protection)
                    for fg_start, fg_end in protection_fg_pairs:
                        if checker.check_fg(fg_start, product_mol) and any(
                            checker.check_fg(fg_end, r) for r in reactant_mols
                        ):
                            print(
                                f"Detected potential deprotection: {fg_end} -> {fg_start} at depth {depth}"
                            )
                            deprotection_events.append(
                                (f"{fg_end} -> {fg_start}", depth)
                            )

                except Exception as e:
                    print(f"Error processing reaction SMILES: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if we have both protection and deprotection events
    if protection_events and deprotection_events:
        # Find the earliest protection and latest deprotection
        min_protection_depth = min(depth for _, depth in protection_events)
        max_deprotection_depth = max(depth for _, depth in deprotection_events)

        # In retrosynthetic traversal, protection should be at a lower depth than deprotection
        # (protection happens earlier in forward synthesis)
        valid_strategy = min_protection_depth < max_deprotection_depth
        print(
            f"Protection depth: {min_protection_depth}, Deprotection depth: {max_deprotection_depth}"
        )
        print(f"Protection-deprotection strategy detected: {valid_strategy}")
        return valid_strategy
    else:
        print(f"Protection-deprotection strategy detected: False (missing events)")
        return False
