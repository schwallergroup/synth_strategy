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
    This function detects if the synthesis involves both protection and deprotection steps.
    """
    protection_found = False
    deprotection_found = False

    # List of protection and deprotection reaction types
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
        "Hydroxyl benzyl deprotection",
        "Carboxyl benzyl deprotection",
        "Cleavage of methoxy ethers to alcohols",
        "Cleavage of alkoxy ethers to alcohols",
        "Ether cleavage to primary alcohol",
        "COOH ethyl deprotection",
        "N-glutarimide deprotection",
        "Phthalimide deprotection",
        "TMS deprotection from alkyne",
    ]

    # Alternative approach: check for functional group changes
    protected_groups = [
        "TMS ether protective group",
        "Silyl protective group",
        "Acetal/Ketal",
        "Boc",
    ]

    unprotected_groups = [
        "Primary alcohol",
        "Secondary alcohol",
        "Tertiary alcohol",
        "Aromatic alcohol",
        "Phenol",
        "Primary amine",
        "Secondary amine",
        "Tertiary amine",
        "Carboxylic acid",
        "Aldehyde",
        "Ketone",
        "Alkyne",
    ]

    def dfs_traverse(node):
        nonlocal protection_found, deprotection_found

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]

                # Check for protection reactions
                for reaction_type in protection_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        print(f"Protection reaction found: {reaction_type}")
                        protection_found = True
                        break

                # Check for deprotection reactions
                for reaction_type in deprotection_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        print(f"Deprotection reaction found: {reaction_type}")
                        deprotection_found = True
                        break

                # If we haven't found a specific reaction type, check for functional group changes
                if not (protection_found and deprotection_found):
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    # Check for protection: unprotected -> protected
                    if not protection_found:
                        for reactant in reactants:
                            # Check if reactant has unprotected group and product has protected group
                            for unprotected in unprotected_groups:
                                if checker.check_fg(unprotected, reactant):
                                    for protected in protected_groups:
                                        if checker.check_fg(
                                            protected, product
                                        ) and not checker.check_fg(protected, reactant):
                                            print(
                                                f"Protection detected: {unprotected} -> {protected}"
                                            )
                                            protection_found = True
                                            break
                                    if protection_found:
                                        break
                            if protection_found:
                                break

                    # Check for deprotection: protected -> unprotected
                    if not deprotection_found:
                        for reactant in reactants:
                            # Check if reactant has protected group and product has unprotected group
                            for protected in protected_groups:
                                if checker.check_fg(protected, reactant):
                                    for unprotected in unprotected_groups:
                                        if checker.check_fg(
                                            unprotected, product
                                        ) and not checker.check_fg(unprotected, reactant):
                                            print(
                                                f"Deprotection detected: {protected} -> {unprotected}"
                                            )
                                            deprotection_found = True
                                            break
                                    if deprotection_found:
                                        break
                            if deprotection_found:
                                break

        # Continue traversing the route
        for child in node.get("children", []):
            dfs_traverse(child)
            # If we've already found both, we can stop traversing
            if protection_found and deprotection_found:
                return

    dfs_traverse(route)
    print(f"Protection found: {protection_found}, Deprotection found: {deprotection_found}")
    return protection_found and deprotection_found
