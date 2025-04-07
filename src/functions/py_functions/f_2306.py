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
    Checks if the synthetic route contains a protection-deprotection sequence.
    """
    # Track protection and deprotection reactions
    protection_reactions = []
    deprotection_reactions = []

    def dfs_traverse(node, depth=0):
        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rxn_smiles = node["metadata"]["rsmi"]

            # Check for protection reactions
            if (
                checker.check_reaction(
                    "Alcohol protection with silyl ethers", rxn_smiles
                )
                or checker.check_reaction("Boc amine protection", rxn_smiles)
                or checker.check_reaction("Boc amine protection explicit", rxn_smiles)
                or checker.check_reaction(
                    "Boc amine protection with Boc anhydride", rxn_smiles
                )
                or checker.check_reaction(
                    "Boc amine protection (ethyl Boc)", rxn_smiles
                )
                or checker.check_reaction(
                    "Boc amine protection of secondary amine", rxn_smiles
                )
                or checker.check_reaction(
                    "Boc amine protection of primary amine", rxn_smiles
                )
                or checker.check_reaction("Protection of carboxylic acid", rxn_smiles)
                or checker.check_reaction("TMS ether protective group", rxn_smiles)
                or checker.check_reaction("Silyl protective group", rxn_smiles)
            ):
                protection_reactions.append((depth, rxn_smiles))
                print(
                    f"Protection reaction found at depth {depth}: {rxn_smiles[:50]}..."
                )

            # Check for deprotection reactions
            if (
                checker.check_reaction(
                    "Alcohol deprotection from silyl ethers", rxn_smiles
                )
                or checker.check_reaction(
                    "Alcohol deprotection from silyl ethers (double)", rxn_smiles
                )
                or checker.check_reaction(
                    "Alcohol deprotection from silyl ethers (diol)", rxn_smiles
                )
                or checker.check_reaction("Boc amine deprotection", rxn_smiles)
                or checker.check_reaction(
                    "Boc amine deprotection of guanidine", rxn_smiles
                )
                or checker.check_reaction(
                    "Boc amine deprotection to NH-NH2", rxn_smiles
                )
                or checker.check_reaction("Deprotection of carboxylic acid", rxn_smiles)
                or checker.check_reaction("TMS deprotection from alkyne", rxn_smiles)
                or checker.check_reaction(
                    "Tert-butyl deprotection of amine", rxn_smiles
                )
                or checker.check_reaction("Hydroxyl benzyl deprotection", rxn_smiles)
                or checker.check_reaction("Carboxyl benzyl deprotection", rxn_smiles)
                or checker.check_reaction(
                    "Cleavage of methoxy ethers to alcohols", rxn_smiles
                )
                or checker.check_reaction(
                    "Cleavage of alkoxy ethers to alcohols", rxn_smiles
                )
                or checker.check_reaction(
                    "Ether cleavage to primary alcohol", rxn_smiles
                )
                or checker.check_reaction("N-glutarimide deprotection", rxn_smiles)
                or checker.check_reaction("Phthalimide deprotection", rxn_smiles)
            ):
                deprotection_reactions.append((depth, rxn_smiles))
                print(
                    f"Deprotection reaction found at depth {depth}: {rxn_smiles[:50]}..."
                )

            # Also check for functional group changes that might indicate protection/deprotection
            reactants = rxn_smiles.split(">")[0].split(".")
            product = rxn_smiles.split(">")[-1]

            # Check for protection patterns
            if any(
                checker.check_fg("Primary alcohol", r) for r in reactants
            ) and checker.check_fg("Ether", product):
                protection_reactions.append((depth, rxn_smiles))
                print(f"Alcohol protection pattern found at depth {depth}")

            if any(
                checker.check_fg("Primary amine", r) for r in reactants
            ) and checker.check_fg("Secondary amide", product):
                protection_reactions.append((depth, rxn_smiles))
                print(f"Amine protection pattern found at depth {depth}")

            # Check for deprotection patterns
            if checker.check_fg("Ether", reactants[0]) and checker.check_fg(
                "Primary alcohol", product
            ):
                deprotection_reactions.append((depth, rxn_smiles))
                print(f"Alcohol deprotection pattern found at depth {depth}")

            if checker.check_fg("Secondary amide", reactants[0]) and checker.check_fg(
                "Primary amine", product
            ):
                deprotection_reactions.append((depth, rxn_smiles))
                print(f"Amine deprotection pattern found at depth {depth}")

        # Recursively traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if we have both protection and deprotection reactions
    if protection_reactions and deprotection_reactions:
        # Check if any protection happens before any deprotection
        valid_sequence = any(
            p_depth < d_depth
            for p_depth, _ in protection_reactions
            for d_depth, _ in deprotection_reactions
        )

        if valid_sequence:
            print(
                f"Protection-deprotection sequence detected: {len(protection_reactions)} protection reactions, {len(deprotection_reactions)} deprotection reactions"
            )
            return True
        else:
            print(f"Protection-deprotection sequence not in correct order")
    else:
        print(
            f"Protection-deprotection sequence not found. Protection reactions: {len(protection_reactions)}, Deprotection reactions: {len(deprotection_reactions)}"
        )

    return False
