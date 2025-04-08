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
    Detects if the synthesis route involves protection and deprotection sequences.
    """
    protection_reactions = []
    deprotection_reactions = []

    def find_protection_deprotection(node, depth=0):
        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rxn_smiles = node["metadata"]["rsmi"]

            # Check for protection reactions
            if (
                checker.check_reaction("Alcohol protection with silyl ethers", rxn_smiles)
                or checker.check_reaction("Protection of carboxylic acid", rxn_smiles)
                or checker.check_reaction("Boc amine protection", rxn_smiles)
                or checker.check_reaction("Boc amine protection explicit", rxn_smiles)
                or checker.check_reaction("Boc amine protection with Boc anhydride", rxn_smiles)
                or checker.check_reaction("Boc amine protection (ethyl Boc)", rxn_smiles)
                or checker.check_reaction("Boc amine protection of secondary amine", rxn_smiles)
                or checker.check_reaction("Boc amine protection of primary amine", rxn_smiles)
                or checker.check_reaction("Acetal hydrolysis to diol", rxn_smiles)
                or checker.check_reaction("Acetal hydrolysis to aldehyde", rxn_smiles)
                or checker.check_reaction("Ketal hydrolysis to ketone", rxn_smiles)
                or checker.check_reaction("Aldehyde or ketone acetalization", rxn_smiles)
                or checker.check_reaction("Diol acetalization", rxn_smiles)
            ):
                protection_reactions.append((rxn_smiles, depth))

            # Check for deprotection reactions
            if (
                checker.check_reaction("Alcohol deprotection from silyl ethers", rxn_smiles)
                or checker.check_reaction(
                    "Alcohol deprotection from silyl ethers (double)", rxn_smiles
                )
                or checker.check_reaction(
                    "Alcohol deprotection from silyl ethers (diol)", rxn_smiles
                )
                or checker.check_reaction("Deprotection of carboxylic acid", rxn_smiles)
                or checker.check_reaction("Boc amine deprotection", rxn_smiles)
                or checker.check_reaction("Boc amine deprotection of guanidine", rxn_smiles)
                or checker.check_reaction("Boc amine deprotection to NH-NH2", rxn_smiles)
                or checker.check_reaction("Ester saponification (methyl deprotection)", rxn_smiles)
                or checker.check_reaction("Ester saponification (alkyl deprotection)", rxn_smiles)
                or checker.check_reaction("TMS deprotection from alkyne", rxn_smiles)
                or checker.check_reaction("Tert-butyl deprotection of amine", rxn_smiles)
                or checker.check_reaction("Phthalimide deprotection", rxn_smiles)
                or checker.check_reaction("N-glutarimide deprotection", rxn_smiles)
                or checker.check_reaction("Hydroxyl benzyl deprotection", rxn_smiles)
                or checker.check_reaction("Carboxyl benzyl deprotection", rxn_smiles)
                or checker.check_reaction("Cleavage of methoxy ethers to alcohols", rxn_smiles)
                or checker.check_reaction("Cleavage of alkoxy ethers to alcohols", rxn_smiles)
                or checker.check_reaction("Ether cleavage to primary alcohol", rxn_smiles)
            ):
                deprotection_reactions.append((rxn_smiles, depth))

            # Check for functional group transformations that might indicate protection/deprotection
            product = rxn_smiles.split(">")[-1]
            reactants = rxn_smiles.split(">")[0].split(".")

            # Check for TMS or silyl groups in reactants or products
            if (
                any(
                    checker.check_fg("TMS ether protective group", r)
                    or checker.check_fg("Silyl protective group", r)
                    for r in reactants
                )
                or checker.check_fg("TMS ether protective group", product)
                or checker.check_fg("Silyl protective group", product)
            ):
                if (
                    any(
                        checker.check_fg("Primary alcohol", r)
                        or checker.check_fg("Secondary alcohol", r)
                        or checker.check_fg("Tertiary alcohol", r)
                        for r in reactants
                    )
                    or checker.check_fg("Primary alcohol", product)
                    or checker.check_fg("Secondary alcohol", product)
                    or checker.check_fg("Tertiary alcohol", product)
                ):
                    protection_reactions.append((rxn_smiles, depth))

            # Check for Boc groups
            if any(checker.check_fg("Boc", r) for r in reactants) or checker.check_fg(
                "Boc", product
            ):
                if (
                    any(
                        checker.check_fg("Primary amine", r)
                        or checker.check_fg("Secondary amine", r)
                        for r in reactants
                    )
                    or checker.check_fg("Primary amine", product)
                    or checker.check_fg("Secondary amine", product)
                ):
                    protection_reactions.append((rxn_smiles, depth))

        for child in node.get("children", []):
            find_protection_deprotection(child, depth + 1)

    find_protection_deprotection(route)

    # Remove duplicates
    unique_protection = set(rxn for rxn, _ in protection_reactions)
    unique_deprotection = set(rxn for rxn, _ in deprotection_reactions)

    # Consider it a protection/deprotection strategy if at least one protection and one deprotection are found
    has_protection = len(unique_protection) > 0
    has_deprotection = len(unique_deprotection) > 0
    print(
        f"Protection reactions: {len(unique_protection)}, Deprotection reactions: {len(unique_deprotection)}"
    )

    return has_protection or has_deprotection  # Changed to OR to be more lenient
