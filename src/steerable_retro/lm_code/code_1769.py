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
    This function detects if the synthetic route involves multiple protection and deprotection steps.
    Returns True if there are at least 2 protection/deprotection steps combined.
    """
    protection_count = 0
    deprotection_count = 0

    # List of common protecting groups
    protecting_groups = [
        "Boc",
        "Silyl protective group",
        "TMS ether protective group",
        "Acetal/Ketal",
        "Carbamic ester",
    ]

    # List of protection/deprotection reaction types
    protection_reactions = [
        "Alcohol protection with silyl ethers",
        "Boc amine protection",
        "Boc amine protection explicit",
        "Boc amine protection with Boc anhydride",
        "Boc amine protection (ethyl Boc)",
        "Boc amine protection of secondary amine",
        "Boc amine protection of primary amine",
        "Protection of carboxylic acid",
        "Aldehyde or ketone acetalization",
        "Diol acetalization",
    ]

    deprotection_reactions = [
        "Alcohol deprotection from silyl ethers",
        "Alcohol deprotection from silyl ethers (double)",
        "Alcohol deprotection from silyl ethers (diol)",
        "Boc amine deprotection",
        "Boc amine deprotection of guanidine",
        "Boc amine deprotection to NH-NH2",
        "Deprotection of carboxylic acid",
        "Acetal hydrolysis to diol",
        "Acetal hydrolysis to aldehyde",
        "Ketal hydrolysis to ketone",
        "Tert-butyl deprotection of amine",
        "TMS deprotection from alkyne",
        "COOH ethyl deprotection",
        "Hydroxyl benzyl deprotection",
        "Carboxyl benzyl deprotection",
        "Cleavage of methoxy ethers to alcohols",
        "Cleavage of alkoxy ethers to alcohols",
        "Ether cleavage to primary alcohol",
        "N-glutarimide deprotection",
        "Phthalimide deprotection",
    ]

    def dfs_traverse(node):
        nonlocal protection_count, deprotection_count

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reaction_found = False

                # Check for protection reactions
                for reaction_type in protection_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        protection_count += 1
                        print(f"Found protection step ({reaction_type}), count: {protection_count}")
                        reaction_found = True
                        break

                # Check for deprotection reactions if no protection reaction was found
                if not reaction_found:
                    for reaction_type in deprotection_reactions:
                        if checker.check_reaction(reaction_type, rsmi):
                            deprotection_count += 1
                            print(
                                f"Found deprotection step ({reaction_type}), count: {deprotection_count}"
                            )
                            reaction_found = True
                            break

                # If no specific reaction type was found, check for protecting group changes
                if not reaction_found:
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    # Check for protection (protecting group in product but not in reactants)
                    product_protecting_groups = []
                    for pg in protecting_groups:
                        if checker.check_fg(pg, product):
                            product_protecting_groups.append(pg)

                    reactant_protecting_groups = []
                    for reactant in reactants:
                        for pg in protecting_groups:
                            if checker.check_fg(pg, reactant):
                                reactant_protecting_groups.append(pg)

                    # New protecting groups added (protection)
                    new_protecting_groups = [
                        pg
                        for pg in product_protecting_groups
                        if pg not in reactant_protecting_groups
                    ]
                    if new_protecting_groups:
                        protection_count += 1
                        print(
                            f"Found protection step (detected {new_protecting_groups[0]}), count: {protection_count}"
                        )

                    # Protecting groups removed (deprotection)
                    removed_protecting_groups = [
                        pg
                        for pg in reactant_protecting_groups
                        if pg not in product_protecting_groups
                    ]
                    if removed_protecting_groups:
                        deprotection_count += 1
                        print(
                            f"Found deprotection step (removed {removed_protecting_groups[0]}), count: {deprotection_count}"
                        )

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    print(f"Total protection steps: {protection_count}")
    print(f"Total deprotection steps: {deprotection_count}")

    # Return True if there are at least 2 protection/deprotection steps combined
    return (protection_count + deprotection_count) >= 2
