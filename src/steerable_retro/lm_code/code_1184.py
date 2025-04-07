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
    This function detects if the synthetic route involves multiple protection steps.
    Returns True if 2 or more protection/deprotection steps are found.
    """
    protection_steps = 0

    def dfs_traverse(node, depth=0):
        nonlocal protection_steps

        # Print node type for debugging
        if node["type"] == "mol":
            print(f"Depth {depth}: Molecule node with SMILES: {node['smiles'][:20]}...")

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            print(f"Depth {depth}: Reaction node with RSMI: {rsmi[:50]}...")

            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # First check for protection/deprotection reactions directly
            # Boc protection
            if (
                checker.check_reaction("Boc amine protection", rsmi)
                or checker.check_reaction("Boc amine protection explicit", rsmi)
                or checker.check_reaction("Boc amine protection with Boc anhydride", rsmi)
                or checker.check_reaction("Boc amine protection (ethyl Boc)", rsmi)
                or checker.check_reaction("Boc amine protection of secondary amine", rsmi)
                or checker.check_reaction("Boc amine protection of primary amine", rsmi)
            ):
                protection_steps += 1
                print(f"Depth {depth}: Boc protection step detected. Total: {protection_steps}")

            # Silyl protection
            elif checker.check_reaction("Alcohol protection with silyl ethers", rsmi):
                protection_steps += 1
                print(f"Depth {depth}: Silyl protection step detected. Total: {protection_steps}")

            # Boc deprotection
            elif (
                checker.check_reaction("Boc amine deprotection", rsmi)
                or checker.check_reaction("Boc amine deprotection of guanidine", rsmi)
                or checker.check_reaction("Boc amine deprotection to NH-NH2", rsmi)
            ):
                protection_steps += 1
                print(f"Depth {depth}: Boc deprotection step detected. Total: {protection_steps}")

            # Silyl deprotection
            elif (
                checker.check_reaction("Alcohol deprotection from silyl ethers", rsmi)
                or checker.check_reaction("Alcohol deprotection from silyl ethers (double)", rsmi)
                or checker.check_reaction("Alcohol deprotection from silyl ethers (diol)", rsmi)
            ):
                protection_steps += 1
                print(f"Depth {depth}: Silyl deprotection step detected. Total: {protection_steps}")

            # Carboxylic acid protection
            elif checker.check_reaction("Protection of carboxylic acid", rsmi):
                protection_steps += 1
                print(
                    f"Depth {depth}: Carboxylic acid protection detected. Total: {protection_steps}"
                )

            # Benzyl deprotection
            elif checker.check_reaction(
                "Hydroxyl benzyl deprotection", rsmi
            ) or checker.check_reaction("Carboxyl benzyl deprotection", rsmi):
                protection_steps += 1
                print(f"Depth {depth}: Benzyl deprotection detected. Total: {protection_steps}")

            # Carboxylic acid deprotection
            elif (
                checker.check_reaction("COOH ethyl deprotection", rsmi)
                or checker.check_reaction("Ester saponification (methyl deprotection)", rsmi)
                or checker.check_reaction("Ester saponification (alkyl deprotection)", rsmi)
                or checker.check_reaction("Deprotection of carboxylic acid", rsmi)
            ):
                protection_steps += 1
                print(
                    f"Depth {depth}: Carboxylic acid deprotection detected. Total: {protection_steps}"
                )

            # TMS alkyne deprotection
            elif checker.check_reaction("TMS deprotection from alkyne", rsmi):
                protection_steps += 1
                print(f"Depth {depth}: TMS alkyne deprotection detected. Total: {protection_steps}")

            # Amine deprotection
            elif checker.check_reaction("Phthalimide deprotection", rsmi) or checker.check_reaction(
                "N-glutarimide deprotection", rsmi
            ):
                protection_steps += 1
                print(f"Depth {depth}: Amine deprotection detected. Total: {protection_steps}")

            # Tert-butyl amine deprotection
            elif checker.check_reaction("Tert-butyl deprotection of amine", rsmi):
                protection_steps += 1
                print(
                    f"Depth {depth}: Tert-butyl amine deprotection detected. Total: {protection_steps}"
                )

            # Acetal/ketal protection
            elif checker.check_reaction(
                "Aldehyde or ketone acetalization", rsmi
            ) or checker.check_reaction("Diol acetalization", rsmi):
                protection_steps += 1
                print(f"Depth {depth}: Acetal/ketal protection detected. Total: {protection_steps}")

            # Acetal/ketal deprotection
            elif (
                checker.check_reaction("Acetal hydrolysis to diol", rsmi)
                or checker.check_reaction("Acetal hydrolysis to aldehyde", rsmi)
                or checker.check_reaction("Ketal hydrolysis to ketone", rsmi)
            ):
                protection_steps += 1
                print(
                    f"Depth {depth}: Acetal/ketal deprotection detected. Total: {protection_steps}"
                )

            # If reaction type checks failed, try functional group checks as fallback
            else:
                # Check for Boc protection by functional groups
                if any(
                    checker.check_fg("Primary amine", r) or checker.check_fg("Secondary amine", r)
                    for r in reactants
                    if r
                ):
                    if checker.check_fg("Boc", product):
                        protection_steps += 1
                        print(
                            f"Depth {depth}: Boc protection detected by FG check. Total: {protection_steps}"
                        )

                # Check for silyl protection by functional groups
                if any(
                    checker.check_fg("Primary alcohol", r)
                    or checker.check_fg("Secondary alcohol", r)
                    or checker.check_fg("Tertiary alcohol", r)
                    for r in reactants
                    if r
                ):
                    if checker.check_fg("Silyl protective group", product) or checker.check_fg(
                        "TMS ether protective group", product
                    ):
                        protection_steps += 1
                        print(
                            f"Depth {depth}: Silyl protection detected by FG check. Total: {protection_steps}"
                        )

                # Check for Boc deprotection by functional groups
                if checker.check_fg("Boc", reactants[0]) if reactants else False:
                    if any(
                        checker.check_fg("Primary amine", p)
                        or checker.check_fg("Secondary amine", p)
                        for p in [product]
                        if p
                    ):
                        protection_steps += 1
                        print(
                            f"Depth {depth}: Boc deprotection detected by FG check. Total: {protection_steps}"
                        )

                # Check for silyl deprotection by functional groups
                if (
                    checker.check_fg("Silyl protective group", reactants[0])
                    or checker.check_fg("TMS ether protective group", reactants[0])
                    if reactants
                    else False
                ):
                    if any(
                        checker.check_fg("Primary alcohol", p)
                        or checker.check_fg("Secondary alcohol", p)
                        or checker.check_fg("Tertiary alcohol", p)
                        for p in [product]
                        if p
                    ):
                        protection_steps += 1
                        print(
                            f"Depth {depth}: Silyl deprotection detected by FG check. Total: {protection_steps}"
                        )

        # Recursively process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Debug output
    print(f"Total protection steps found: {protection_steps}")
    result = protection_steps >= 2
    print(f"Returning: {result}")

    return result
