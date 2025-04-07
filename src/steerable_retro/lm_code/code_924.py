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
    This function detects if the synthesis uses TMS protection for a terminal alkyne.
    """
    has_tms_protection = False

    def dfs_traverse(node):
        nonlocal has_tms_protection

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            print(f"Analyzing reaction: {rsmi}")

            # Check for TMS protection/deprotection reactions
            is_silyl_protection = checker.check_reaction(
                "Alcohol protection with silyl ethers", rsmi
            )
            is_silyl_deprotection = checker.check_reaction(
                "Alcohol deprotection from silyl ethers", rsmi
            ) or checker.check_reaction("TMS deprotection from alkyne", rsmi)

            # Forward direction: Terminal alkyne → TMS-protected alkyne
            if any(checker.check_fg("Alkyne", r) for r in reactants_smiles) and checker.check_fg(
                "Alkyne", product_smiles
            ):
                # Protection: reactant has terminal alkyne, product has TMS-protected alkyne
                reactant_has_terminal_alkyne = any(
                    checker.check_fg("Alkyne", r)
                    and not checker.check_fg("TMS ether protective group", r)
                    for r in reactants_smiles
                )
                product_has_tms_alkyne = checker.check_fg(
                    "TMS ether protective group", product_smiles
                ) and checker.check_fg("Alkyne", product_smiles)

                if reactant_has_terminal_alkyne and product_has_tms_alkyne:
                    print(f"Found TMS protection of alkyne (forward): {rsmi}")
                    has_tms_protection = True

                # Deprotection: reactant has TMS-protected alkyne, product has terminal alkyne
                reactant_has_tms_alkyne = any(
                    checker.check_fg("TMS ether protective group", r)
                    and checker.check_fg("Alkyne", r)
                    for r in reactants_smiles
                )
                product_has_terminal_alkyne = checker.check_fg(
                    "Alkyne", product_smiles
                ) and not checker.check_fg("TMS ether protective group", product_smiles)

                if reactant_has_tms_alkyne and product_has_terminal_alkyne:
                    print(f"Found TMS deprotection of alkyne (forward): {rsmi}")
                    has_tms_protection = True

            # Additional check for TMS groups that might not be recognized as protective groups
            if "Si" in rsmi and any(
                checker.check_fg("Alkyne", r) for r in reactants_smiles + [product_smiles]
            ):
                # Look for trimethylsilyl group pattern in SMILES
                tms_pattern = "C[Si](C)(C)"

                # Protection: TMS appears in product but not in reactants
                if tms_pattern not in "".join(reactants_smiles) and tms_pattern in product_smiles:
                    if any(checker.check_fg("Alkyne", r) for r in reactants_smiles):
                        print(f"Found TMS protection via Si pattern: {rsmi}")
                        has_tms_protection = True

                # Deprotection: TMS appears in reactants but not in product
                if (
                    any(tms_pattern in r for r in reactants_smiles)
                    and tms_pattern not in product_smiles
                ):
                    if checker.check_fg("Alkyne", product_smiles):
                        print(f"Found TMS deprotection via Si pattern: {rsmi}")
                        has_tms_protection = True

        # Continue DFS traversal
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    print(f"TMS protection strategy detected: {has_tms_protection}")
    return has_tms_protection
