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
    This function detects if the synthesis uses epoxide opening to form a tertiary alcohol.
    """
    epoxide_opening_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal epoxide_opening_detected

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if this is an epoxide opening reaction
                if checker.check_reaction("Ring opening of epoxide with amine", rsmi):
                    print(f"Epoxide opening reaction detected at depth {depth}")

                    # Check if product has a tertiary alcohol
                    if checker.check_fg("Tertiary alcohol", product):
                        print(f"Tertiary alcohol found in product at depth {depth}")
                        epoxide_opening_detected = True
                        return

                # If not a standard epoxide opening reaction, check manually
                for reactant in reactants:
                    if checker.check_ring("oxirane", reactant):
                        print(f"Epoxide found in reactant at depth {depth}")

                        # Check if product has tertiary alcohol
                        if checker.check_fg("Tertiary alcohol", product):
                            print(f"Tertiary alcohol found in product at depth {depth}")

                            # If we found an epoxide in reactant and tertiary alcohol in product,
                            # this is sufficient evidence of epoxide opening
                            print(f"Epoxide opening to tertiary alcohol detected at depth {depth}")
                            epoxide_opening_detected = True
                            return

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return epoxide_opening_detected
