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
    Detects the use of Boc (tert-butoxycarbonyl) protection strategy for amines.
    """
    has_boc_protection = False

    def dfs_traverse(node, depth=0):
        nonlocal has_boc_protection

        # Check for molecules that already contain Boc-protected amines
        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            if checker.check_fg("Carbamic ester", mol_smiles):
                has_boc_protection = True
                print(
                    f"Detected molecule with Boc-protected amine (carbamic ester) at depth {depth}"
                )

        # Check for Boc protection/deprotection reactions
        elif node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for Boc protection reactions
                if (
                    checker.check_reaction("Boc amine protection", rsmi)
                    or checker.check_reaction("Boc amine protection explicit", rsmi)
                    or checker.check_reaction("Boc amine protection with Boc anhydride", rsmi)
                    or checker.check_reaction("Boc amine protection (ethyl Boc)", rsmi)
                    or checker.check_reaction("Boc amine protection of secondary amine", rsmi)
                    or checker.check_reaction("Boc amine protection of primary amine", rsmi)
                ):
                    has_boc_protection = True
                    print(f"Detected Boc protection reaction at depth {depth}")

                # Check for Boc deprotection reactions
                elif (
                    checker.check_reaction("Boc amine deprotection", rsmi)
                    or checker.check_reaction("Boc amine deprotection of guanidine", rsmi)
                    or checker.check_reaction("Boc amine deprotection to NH-NH2", rsmi)
                    or checker.check_reaction("Tert-butyl deprotection of amine", rsmi)
                ):
                    has_boc_protection = True
                    print(f"Detected Boc deprotection reaction at depth {depth}")

                # Check for Boc protection reactions by examining functional group changes
                else:
                    for r in reactants:
                        if checker.check_fg("Primary amine", r) or checker.check_fg(
                            "Secondary amine", r
                        ):
                            # Check if product contains a carbamic ester (characteristic of Boc protection)
                            if checker.check_fg("Carbamic ester", product):
                                has_boc_protection = True
                                print(f"Detected Boc-protected amine formation at depth {depth}")
                                break
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Boc protection strategy: {has_boc_protection}")
    return has_boc_protection
