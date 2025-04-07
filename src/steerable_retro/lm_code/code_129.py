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
    Detects if the synthesis route includes a carboxylic acid protection-deprotection
    cycle using benzyl groups.
    """
    # Track if we've found protection and deprotection steps
    found_protection = False
    found_deprotection = False

    def dfs(node, depth=0):
        nonlocal found_protection, found_deprotection

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]

            # Check for benzyl protection of carboxylic acid
            if checker.check_reaction("Protection of carboxylic acid", rsmi):
                # Verify it's specifically benzyl protection
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if reactant has carboxylic acid and product has benzyl ester
                for reactant in reactants:
                    if checker.check_fg("Carboxylic acid", reactant):
                        # Check for benzyl group in product
                        if checker.check_fg("Ester", product) and "c1ccccc1" in product:
                            print(f"Found benzyl protection of carboxylic acid at depth {depth}")
                            found_protection = True

            # Check for benzyl deprotection
            if checker.check_reaction("Carboxyl benzyl deprotection", rsmi):
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Verify it's a benzyl deprotection to carboxylic acid
                for reactant in reactants:
                    if "c1ccccc1" in reactant and checker.check_fg("Ester", reactant):
                        if checker.check_fg("Carboxylic acid", product):
                            print(f"Found benzyl deprotection to carboxylic acid at depth {depth}")
                            found_deprotection = True

        # Continue DFS traversal
        for child in node.get("children", []):
            dfs(child, depth + 1)

    # Start DFS traversal
    dfs(route)

    return found_protection and found_deprotection
