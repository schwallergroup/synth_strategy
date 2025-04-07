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
    Detects if the synthesis route involves Boc protection followed by deprotection.
    """
    boc_protection_found = False
    boc_deprotection_found = False

    def dfs_traverse(node, depth=0):
        nonlocal boc_protection_found, boc_deprotection_found

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
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
                    print(f"Found Boc protection reaction at depth {depth}: {rsmi}")
                    boc_protection_found = True

                # Check for Boc deprotection reactions
                if (
                    checker.check_reaction("Boc amine deprotection", rsmi)
                    or checker.check_reaction("Boc amine deprotection of guanidine", rsmi)
                    or checker.check_reaction("Boc amine deprotection to NH-NH2", rsmi)
                    or checker.check_reaction("Tert-butyl deprotection of amine", rsmi)
                ):
                    print(f"Found Boc deprotection reaction at depth {depth}: {rsmi}")
                    boc_deprotection_found = True

                # Additional check: verify Boc group is actually being added/removed
                if not boc_protection_found:
                    # Check if any reactant doesn't have Boc but product does
                    reactant_has_boc = all(
                        checker.check_fg("Boc", r) for r in reactants if r.strip()
                    )
                    product_has_boc = checker.check_fg("Boc", product)
                    if not reactant_has_boc and product_has_boc:
                        print(f"Found Boc addition (FG check) at depth {depth}: {rsmi}")
                        boc_protection_found = True

                if not boc_deprotection_found:
                    # Check if any reactant has Boc but product doesn't
                    reactant_has_boc = any(
                        checker.check_fg("Boc", r) for r in reactants if r.strip()
                    )
                    product_has_boc = checker.check_fg("Boc", product)
                    if reactant_has_boc and not product_has_boc:
                        print(f"Found Boc removal (FG check) at depth {depth}: {rsmi}")
                        boc_deprotection_found = True

        # Recursively process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    print(f"Protection found: {boc_protection_found}, Deprotection found: {boc_deprotection_found}")
    return boc_protection_found and boc_deprotection_found
