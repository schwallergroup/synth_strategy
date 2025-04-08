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
    This function detects if the synthesis uses a ketal protection/deprotection strategy
    for a carbonyl group.
    """
    found_ketal_formation = False
    found_ketal_deprotection = False

    def dfs_traverse(node, depth=0):
        nonlocal found_ketal_formation, found_ketal_deprotection

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for ketal formation (protection)
            if checker.check_reaction("Aldehyde or ketone acetalization", rsmi):
                # Verify ketone in reactants and ketal in product
                if any(checker.check_fg("Ketone", r) for r in reactants) and checker.check_fg(
                    "Acetal/Ketal", product
                ):
                    print(f"Found ketal formation at depth {depth}, rsmi: {rsmi}")
                    found_ketal_formation = True

            # Check for ketal deprotection
            if checker.check_reaction("Ketal hydrolysis to ketone", rsmi):
                # Verify ketal in reactants and ketone in product
                if any(checker.check_fg("Acetal/Ketal", r) for r in reactants) and checker.check_fg(
                    "Ketone", product
                ):
                    print(f"Found ketal deprotection at depth {depth}, rsmi: {rsmi}")
                    found_ketal_deprotection = True

            # Fallback detection if reaction checkers don't identify the reactions
            if not found_ketal_formation:
                # Check if reactants contain ketone and product contains ketal
                if any(checker.check_fg("Ketone", r) for r in reactants) and checker.check_fg(
                    "Acetal/Ketal", product
                ):
                    print(
                        f"Found potential ketal formation at depth {depth} (fallback method), rsmi: {rsmi}"
                    )
                    found_ketal_formation = True

            if not found_ketal_deprotection:
                # Check if reactants contain ketal and product contains ketone
                if any(checker.check_fg("Acetal/Ketal", r) for r in reactants) and checker.check_fg(
                    "Ketone", product
                ):
                    print(
                        f"Found potential ketal deprotection at depth {depth} (fallback method), rsmi: {rsmi}"
                    )
                    found_ketal_deprotection = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)

    # Return True only if both protection and deprotection are found
    result = found_ketal_formation and found_ketal_deprotection
    print(f"Ketal protection strategy detected: {result}")
    print(f"Found ketal formation: {found_ketal_formation}")
    print(f"Found ketal deprotection: {found_ketal_deprotection}")

    return result
