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
    This function detects a strategy where an alkyne intermediate is used
    to form a heterocycle (specifically looking for terminal alkyne → heterocycle).
    """
    # Define heterocycle ring types to check
    heterocycle_rings = [
        "pyrrole",
        "pyridine",
        "pyrazole",
        "imidazole",
        "oxazole",
        "thiazole",
        "triazole",
        "tetrazole",
        "indole",
        "benzimidazole",
        "benzoxazole",
        "benzothiazole",
        "indazole",
        "isoxazole",
        "isothiazole",
        "oxadiazole",
        "thiadiazole",
        "furan",
        "thiophene",
        "pyran",
        "dioxane",
    ]

    # Track reactions in the synthetic route
    reaction_sequence = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Store reaction information with depth
            reaction_info = {
                "rsmi": rsmi,
                "reactants": reactants,
                "product": product,
                "depth": depth,
                "has_alkyne_reactant": False,
                "forms_heterocycle": False,
                "heterocycle_type": None,
            }

            try:
                # Check for terminal alkyne in reactants
                for reactant in reactants:
                    if checker.check_fg("Alkyne", reactant):
                        reaction_info["has_alkyne_reactant"] = True
                        print(f"Terminal alkyne detected in reaction at depth {depth}")
                        break

                # Check for heterocycle formation in product
                product_has_heterocycle = False
                product_heterocycle_type = None

                # Check if product contains a heterocycle that wasn't in the reactants
                for ring_type in heterocycle_rings:
                    if checker.check_ring(ring_type, product):
                        # Check if any reactant already had this heterocycle
                        reactants_have_heterocycle = any(
                            checker.check_ring(ring_type, reactant) for reactant in reactants
                        )

                        if not reactants_have_heterocycle:
                            product_has_heterocycle = True
                            product_heterocycle_type = ring_type
                            print(f"Heterocycle formation ({ring_type}) detected at depth {depth}")
                            break

                reaction_info["forms_heterocycle"] = product_has_heterocycle
                reaction_info["heterocycle_type"] = product_heterocycle_type

                # Add to our sequence
                reaction_sequence.append(reaction_info)

            except Exception as e:
                print(f"Error analyzing reaction: {e}")

        # Traverse children with increased depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Analyze the reaction sequence to find alkyne → heterocycle pattern
    # Sort by depth to get chronological order (higher depth = earlier in synthesis)
    reaction_sequence.sort(key=lambda x: x["depth"], reverse=True)

    # Look for the pattern: a reaction with alkyne followed by a reaction forming heterocycle
    found_alkyne_step = False
    alkyne_depth = -1

    for reaction in reaction_sequence:
        if reaction["has_alkyne_reactant"] and not found_alkyne_step:
            found_alkyne_step = True
            alkyne_depth = reaction["depth"]
            print(f"Found alkyne intermediate at depth {alkyne_depth}")

        # Check if we found a heterocycle formation after the alkyne step
        if found_alkyne_step and reaction["forms_heterocycle"] and reaction["depth"] < alkyne_depth:
            print(
                f"Found heterocycle formation at depth {reaction['depth']} after alkyne at depth {alkyne_depth}"
            )
            print(f"Heterocycle type: {reaction['heterocycle_type']}")
            return True

    return False
