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
    This function detects a late-stage esterification strategy where an acid chloride
    and an alcohol are combined to form an ester in the final synthetic step.
    """
    print("Starting late_stage_esterification_strategy analysis")
    esterification_found = False

    def dfs_traverse(node, depth=0):
        nonlocal esterification_found

        print(f"Traversing node of type {node['type']} at depth {depth}")

        # Check if this is a reaction node
        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            # Only check the final, penultimate, or antepenultimate reaction (depth 0, 1, or 2)
            if depth <= 2:
                print(f"Checking reaction at depth {depth}")
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Reaction SMILES: {rsmi}")

                # Check for various esterification reactions
                is_esterification = False

                # Check specific esterification reaction types
                if checker.check_reaction("Schotten-Baumann to ester", rsmi):
                    print("Schotten-Baumann esterification reaction detected")
                    is_esterification = True
                elif checker.check_reaction("Esterification of Carboxylic Acids", rsmi):
                    print("Carboxylic acid esterification reaction detected")
                    is_esterification = True
                elif checker.check_reaction("Transesterification", rsmi):
                    print("Transesterification reaction detected")
                    is_esterification = True
                elif checker.check_reaction("Oxidative esterification of primary alcohols", rsmi):
                    print("Oxidative esterification reaction detected")
                    is_esterification = True
                elif checker.check_reaction(
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_OS",
                    rsmi,
                ):
                    print("Acylation of O-nucleophiles reaction detected")
                    is_esterification = True

                # If no specific reaction type detected, check for general pattern
                if not is_esterification:
                    # Check if product contains an ester that wasn't in the reactants
                    if checker.check_fg("Ester", product):
                        print("Ester found in product, checking if it's a new formation")
                        ester_in_reactants = any(checker.check_fg("Ester", r) for r in reactants)
                        if not ester_in_reactants:
                            print("New ester formation detected")
                            is_esterification = True

                if is_esterification:
                    # Verify reactants contain acyl halide/carboxylic acid and alcohol
                    acyl_source_found = False
                    alcohol_found = False

                    for reactant in reactants:
                        try:
                            # Check for acyl sources
                            if checker.check_fg("Acyl halide", reactant):
                                acyl_source_found = True
                                print(f"Acyl halide found in reactant: {reactant}")
                            elif checker.check_fg("Carboxylic acid", reactant):
                                acyl_source_found = True
                                print(f"Carboxylic acid found in reactant: {reactant}")
                            elif checker.check_fg("Anhydride", reactant):
                                acyl_source_found = True
                                print(f"Anhydride found in reactant: {reactant}")

                            # Check for alcohols - include all types and a simple OH check
                            if (
                                checker.check_fg("Primary alcohol", reactant)
                                or checker.check_fg("Secondary alcohol", reactant)
                                or checker.check_fg("Tertiary alcohol", reactant)
                                or checker.check_fg("Aromatic alcohol", reactant)
                                or checker.check_fg("Phenol", reactant)
                                or "[OH]" in reactant
                                or "OH" in reactant
                            ):
                                alcohol_found = True
                                print(f"Alcohol found in reactant: {reactant}")
                        except Exception as e:
                            print(f"Error checking reactant {reactant}: {e}")

                    # Verify product contains ester
                    try:
                        if checker.check_fg("Ester", product):
                            print(f"Ester found in product: {product}")

                            if acyl_source_found and alcohol_found:
                                esterification_found = True
                                print(
                                    f"Late-stage esterification strategy confirmed at depth {depth}"
                                )
                        else:
                            print(f"No ester found in product: {product}")
                    except Exception as e:
                        print(f"Error checking product {product}: {e}")

        # Continue traversing with incremented depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    print(f"Esterification strategy found: {esterification_found}")
    return esterification_found
