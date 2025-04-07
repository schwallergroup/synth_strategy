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
    Detects if the synthesis involves late-stage alkyne formation followed by
    aryl-alkyne disconnection (Sonogashira-type pattern).
    """
    alkyne_formation_found = False
    aryl_alkyne_disconnection_found = False

    def dfs_traverse(node, depth=0):
        nonlocal alkyne_formation_found, aryl_alkyne_disconnection_found

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Examining reaction at depth {depth}: {rsmi}")

                # Check for Sonogashira reaction at any depth
                sonogashira_variants = [
                    "Sonogashira acetylene_aryl halide",
                    "Sonogashira alkyne_aryl halide",
                    "Sonogashira acetylene_aryl OTf",
                    "Sonogashira alkyne_aryl OTf",
                    "Sonogashira acetylene_alkenyl halide",
                    "Sonogashira alkyne_alkenyl halide",
                    "Sonogashira acetylene_alkenyl OTf",
                    "Sonogashira alkyne_alkenyl OTf",
                    "Sonogashira acetylene_acyl halide",
                    "Sonogashira alkyne_acyl halide",
                ]

                for variant in sonogashira_variants:
                    if checker.check_reaction(variant, rsmi):
                        aryl_alkyne_disconnection_found = True
                        print(f"Found {variant} at depth {depth}")
                        break

                # If no Sonogashira reaction was detected, check manually for the pattern
                if not aryl_alkyne_disconnection_found:
                    if checker.check_fg("Alkyne", product):
                        # Check if one reactant has aryl-halide and another has terminal alkyne
                        aryl_halide_found = False
                        terminal_alkyne_found = False

                        for reactant in reactants:
                            if checker.check_fg("Aromatic halide", reactant):
                                aryl_halide_found = True
                                print(f"Reactant contains aromatic halide: {reactant}")

                            if checker.check_fg("Alkyne", reactant):
                                terminal_alkyne_found = True
                                print(f"Reactant contains alkyne: {reactant}")

                        if aryl_halide_found and terminal_alkyne_found:
                            aryl_alkyne_disconnection_found = True
                            print(f"Found aryl-alkyne disconnection pattern at depth {depth}")

                # Check for alkyne formation (product has alkyne but reactants don't)
                if checker.check_fg("Alkyne", product):
                    print(f"Product contains alkyne: {product}")

                    # Check if reactants don't have alkyne
                    reactants_have_alkyne = False
                    for reactant in reactants:
                        if checker.check_fg("Alkyne", reactant):
                            reactants_have_alkyne = True
                            print(f"Reactant contains alkyne: {reactant}")
                            break

                    if not reactants_have_alkyne:
                        alkyne_formation_found = True
                        print(f"Found alkyne formation at depth {depth}")

            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"alkyne_formation_found: {alkyne_formation_found}")
    print(f"aryl_alkyne_disconnection_found: {aryl_alkyne_disconnection_found}")

    # Return True if Sonogashira-type disconnection is found
    # In retrosynthetic analysis, we're looking for the Sonogashira pattern
    return aryl_alkyne_disconnection_found
