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
    This function detects a synthetic strategy involving thiazole ring formation
    via a thioamide intermediate.
    """
    # Track if we found the key intermediates and transformations
    found_thioamide = False
    found_thiazole_formation = False

    def dfs_traverse(node):
        nonlocal found_thioamide, found_thiazole_formation

        if node["type"] == "mol":
            # Check if this molecule contains a thioamide group
            if checker.check_fg("Thioamide", node["smiles"]):
                found_thioamide = True
                print(f"Found thioamide intermediate: {node['smiles']}")

        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for thioamide in reactants
            for reactant in reactants:
                if checker.check_fg("Thioamide", reactant):
                    found_thioamide = True
                    print(f"Found thioamide in reaction: {reactant}")

            # Check for thiazole formation
            if checker.check_ring("thiazole", product):
                # Check if any reactant has thioamide but not thiazole
                for reactant in reactants:
                    if checker.check_fg("Thioamide", reactant) and not checker.check_ring(
                        "thiazole", reactant
                    ):
                        found_thiazole_formation = True
                        print(f"Found thiazole formation from thioamide: {rsmi}")
                        break

            # Alternative: check if this is a thiazole formation reaction
            if checker.check_reaction("{thiazole}", rsmi):
                print(f"Found thiazole formation reaction: {rsmi}")
                # Verify that a thioamide is converted to a thiazole
                for reactant in reactants:
                    if checker.check_fg("Thioamide", reactant) and not checker.check_ring(
                        "thiazole", reactant
                    ):
                        if checker.check_ring("thiazole", product):
                            found_thiazole_formation = True
                            print(f"Confirmed thiazole formation from thioamide: {rsmi}")
                            break

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if both conditions are met
    return found_thioamide and found_thiazole_formation
