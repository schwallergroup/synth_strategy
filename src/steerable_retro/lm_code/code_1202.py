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
    Detects synthesis routes that form benzothiazole rings from thiourea intermediates.
    """
    # Track if we found the strategy
    found_thiourea = False
    found_benzothiazole_formation = False

    def dfs_traverse(node):
        nonlocal found_thiourea, found_benzothiazole_formation

        if node["type"] == "mol":
            # Check for thiourea in molecules
            if checker.check_fg("Thiourea", node["smiles"]):
                found_thiourea = True
                print(f"Found thiourea intermediate: {node['smiles']}")

        elif node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"].get("rsmi", "")
            if not rsmi:
                return

            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check if this is a benzothiazole formation reaction
            if checker.check_reaction("benzothiazole", rsmi) or checker.check_reaction(
                "{benzothiazole}", rsmi
            ):
                print(f"Found benzothiazole formation reaction: {rsmi}")

                # Check if product contains benzothiazole ring
                if checker.check_ring("benzothiazole", product_smiles):
                    print(f"Product contains benzothiazole ring: {product_smiles}")

                    # Check if any reactant has thiourea
                    for reactant in reactants_smiles:
                        if checker.check_fg("Thiourea", reactant):
                            found_benzothiazole_formation = True
                            print(f"Found benzothiazole formation from thiourea")
                            break

            # Alternative check: look for benzothiazole formation without specific reaction type
            if not found_benzothiazole_formation:
                # Check if product contains benzothiazole ring
                if checker.check_ring("benzothiazole", product_smiles):
                    # Check if any reactant has thiourea but not benzothiazole
                    for reactant in reactants_smiles:
                        if checker.check_fg("Thiourea", reactant) and not checker.check_ring(
                            "benzothiazole", reactant
                        ):
                            found_benzothiazole_formation = True
                            print(f"Found benzothiazole formation from thiourea (pattern-based)")
                            break

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if both conditions are met
    return found_thiourea and found_benzothiazole_formation
