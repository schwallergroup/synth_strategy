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
    This function detects a synthetic strategy involving TBDMS protection and deprotection
    of alcohols during the synthesis.
    """
    has_tbdms_protection = False
    has_tbdms_deprotection = False

    def dfs_traverse(node, depth=0):
        nonlocal has_tbdms_protection, has_tbdms_deprotection

        # Process reaction nodes
        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            print(f"Examining reaction at depth {depth}: {rsmi}")

            # Check for TBDMS protection
            reactants = rsmi.split(">")[0]
            product = rsmi.split(">")[-1]

            # Check for silyl protection reaction
            if checker.check_reaction("Alcohol protection with silyl ethers", rsmi):
                has_tbdms_protection = True
                print(f"Found TBDMS protection reaction at depth {depth}: {rsmi}")
            # Alternative check: silyl chloride + alcohol → silyl ether
            elif (
                "Si" in reactants
                and "Cl" in reactants
                and "OH" in reactants
                and "Si" in product
                and "O" in product
            ):
                # Check if reactants contain silyl chloride and alcohol
                for reactant in reactants.split("."):
                    if (
                        checker.check_fg("Silyl protective group", reactant)
                        or "Si" in reactant
                        and "Cl" in reactant
                    ):
                        # Check if product contains silyl ether
                        if checker.check_fg("Silyl protective group", product) or checker.check_fg(
                            "TMS ether protective group", product
                        ):
                            has_tbdms_protection = True
                            print(f"Found silyl protection reaction at depth {depth}: {rsmi}")
                            break

            # Check for TBDMS deprotection
            if (
                checker.check_reaction("Alcohol deprotection from silyl ethers", rsmi)
                or checker.check_reaction("Alcohol deprotection from silyl ethers (double)", rsmi)
                or checker.check_reaction("Alcohol deprotection from silyl ethers (diol)", rsmi)
            ):
                has_tbdms_deprotection = True
                print(f"Found TBDMS deprotection reaction at depth {depth}: {rsmi}")
            # Alternative check: silyl ether → alcohol
            elif "Si" in reactants and "O" in reactants and "OH" in product:
                # Check if reactants contain silyl ether
                for reactant in reactants.split("."):
                    if checker.check_fg("Silyl protective group", reactant) or checker.check_fg(
                        "TMS ether protective group", reactant
                    ):
                        has_tbdms_deprotection = True
                        print(f"Found silyl deprotection reaction at depth {depth}: {rsmi}")
                        break

        # Process molecule nodes (optional additional check)
        elif node["type"] == "mol" and "smiles" in node:
            mol_smiles = node["smiles"]
            if checker.check_fg("TMS ether protective group", mol_smiles) or checker.check_fg(
                "Silyl protective group", mol_smiles
            ):
                print(f"Found molecule with silyl ether at depth {depth}: {mol_smiles}")

        # Recursively traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    print("Starting traversal of synthesis route")
    dfs_traverse(route)

    # Report findings
    print(f"Protection found: {has_tbdms_protection}, Deprotection found: {has_tbdms_deprotection}")

    # Return True if both protection and deprotection are found
    return has_tbdms_protection and has_tbdms_deprotection
