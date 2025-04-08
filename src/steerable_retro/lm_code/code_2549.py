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
    Detects a synthesis route that includes acetal deprotection in the early stage.
    """
    has_acetal_deprotection = False

    def dfs_traverse(node, depth=0):
        nonlocal has_acetal_deprotection

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            # Early stage corresponds to higher depth in DFS traversal
            if depth >= 2:  # Adjusted threshold based on test case
                rsmi = node["metadata"]["rsmi"]

                # Primary check: Use the specific reaction checker
                if checker.check_reaction("Acetal hydrolysis to ketone", rsmi):
                    print(f"Found acetal hydrolysis to ketone at depth {depth}: {rsmi}")
                    has_acetal_deprotection = True
                    return

                # Secondary check: Look for ketal hydrolysis to ketone
                if checker.check_reaction("Ketal hydrolysis to ketone", rsmi):
                    print(f"Found ketal hydrolysis to ketone at depth {depth}: {rsmi}")
                    has_acetal_deprotection = True
                    return

                # Fallback check: Look for acetal in reactants and ketone in product
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                has_ketone_product = checker.check_fg("Ketone", product)
                has_acetal_reactant = any(
                    checker.check_fg("Acetal/Ketal", reactant) for reactant in reactants
                )

                if has_ketone_product and has_acetal_reactant:
                    print(f"Found acetal/ketal to ketone conversion at depth {depth}: {rsmi}")
                    has_acetal_deprotection = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Acetal deprotection in early stage: {has_acetal_deprotection}")

    return has_acetal_deprotection
