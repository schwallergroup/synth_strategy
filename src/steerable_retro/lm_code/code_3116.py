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
    This function detects if the synthesis route includes hydration of an alkene
    to form an alcohol.
    """
    hydration_found = False

    def dfs_traverse(node):
        nonlocal hydration_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this is a known alkene hydration reaction
            if checker.check_reaction(
                "Markovnikov alkene hydration to alcohol", rsmi
            ) or checker.check_reaction("anti-Markovnikov alkene hydration to alcohol", rsmi):
                print(f"Alkene hydration reaction detected: {rsmi}")
                hydration_found = True
            else:
                # Fallback: Check for alkene in reactants and alcohol in products
                reactant_has_alkene = False
                product_has_alcohol = False

                for reactant in reactants:
                    if checker.check_fg("Alkene", reactant) or checker.check_fg("Vinyl", reactant):
                        print(f"Alkene found in reactant: {reactant}")
                        reactant_has_alkene = True
                        break

                if reactant_has_alkene:
                    if (
                        checker.check_fg("Primary alcohol", product)
                        or checker.check_fg("Secondary alcohol", product)
                        or checker.check_fg("Tertiary alcohol", product)
                    ):
                        print(f"Alcohol found in product: {product}")
                        product_has_alcohol = True

                if reactant_has_alkene and product_has_alcohol:
                    # Additional check to ensure it's a hydration reaction
                    # Look for common hydration reagents or conditions in the reaction
                    if "H2O" in rsmi or "OH" in rsmi:
                        print(f"Alkene hydration detected through pattern matching: {rsmi}")
                        hydration_found = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(f"Hydration strategy found: {hydration_found}")
    return hydration_found
