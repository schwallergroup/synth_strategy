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
    Detects if the synthetic route employs hydroxyl protection as a key strategy.
    """
    hydroxyl_protections = 0

    def dfs_traverse(node):
        nonlocal hydroxyl_protections

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for hydroxyl protection reactions
                if any(
                    checker.check_fg("Primary alcohol", r)
                    or checker.check_fg("Secondary alcohol", r)
                    or checker.check_fg("Tertiary alcohol", r)
                    or checker.check_fg("Phenol", r)
                    or checker.check_fg("Aromatic alcohol", r)
                    for r in reactants
                ):

                    # Check for common protection reactions
                    if (
                        checker.check_reaction("Alcohol protection with silyl ethers", rsmi)
                        or checker.check_fg("TMS ether protective group", product)
                        or checker.check_fg("Silyl protective group", product)
                        or (
                            checker.check_reaction("Aldehyde or ketone acetalization", rsmi)
                            and any(
                                checker.check_fg("Primary alcohol", r)
                                or checker.check_fg("Secondary alcohol", r)
                                or checker.check_fg("Tertiary alcohol", r)
                                for r in reactants
                            )
                        )
                        or checker.check_reaction("Diol acetalization", rsmi)
                        or (
                            checker.check_reaction("Williamson Ether Synthesis", rsmi)
                            and not checker.check_reaction(
                                "Williamson Ether Synthesis (intra to epoxy)", rsmi
                            )
                        )
                    ):

                        hydroxyl_protections += 1
                        print(f"Found hydroxyl protection: {rsmi}")

                # Check for hydroxyl deprotection reactions
                if (
                    checker.check_reaction("Alcohol deprotection from silyl ethers", rsmi)
                    or checker.check_reaction(
                        "Alcohol deprotection from silyl ethers (double)", rsmi
                    )
                    or checker.check_reaction("Alcohol deprotection from silyl ethers (diol)", rsmi)
                    or checker.check_reaction("Acetal hydrolysis to diol", rsmi)
                    or checker.check_reaction("Acetal hydrolysis to aldehyde", rsmi)
                    or checker.check_reaction("Ketal hydrolysis to ketone", rsmi)
                    or checker.check_reaction("Cleavage of methoxy ethers to alcohols", rsmi)
                    or checker.check_reaction("Cleavage of alkoxy ethers to alcohols", rsmi)
                    or checker.check_reaction("Ether cleavage to primary alcohol", rsmi)
                    or checker.check_reaction("Hydroxyl benzyl deprotection", rsmi)
                ):

                    # Verify that a hydroxyl group is produced
                    if any(
                        checker.check_fg("Primary alcohol", product)
                        or checker.check_fg("Secondary alcohol", product)
                        or checker.check_fg("Tertiary alcohol", product)
                        or checker.check_fg("Phenol", product)
                        or checker.check_fg("Aromatic alcohol", product)
                    ):

                        hydroxyl_protections += 1
                        print(f"Found hydroxyl deprotection: {rsmi}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(f"Total hydroxyl protections: {hydroxyl_protections}")
    return hydroxyl_protections >= 1
