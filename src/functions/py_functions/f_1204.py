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
    Detects synthesis routes that use thiourea intermediates as precursors for heterocycle formation.
    """
    # Track if we found the strategy
    found_thiourea_intermediate = False
    thiourea_used_for_heterocycle = False

    def dfs_traverse(node):
        nonlocal found_thiourea_intermediate, thiourea_used_for_heterocycle

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"].get("rsmi", "")
            if not rsmi:
                return

            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check if any reactant has thiourea
            has_thiourea_in_reactants = False
            thiourea_reactant = None

            for reactant in reactants_smiles:
                try:
                    if checker.check_fg("Thiourea", reactant):
                        has_thiourea_in_reactants = True
                        found_thiourea_intermediate = True
                        thiourea_reactant = reactant
                        print(f"Found thiourea intermediate: {reactant}")
                        break
                except Exception as e:
                    print(f"Error checking thiourea in reactant: {e}")
                    continue

            # If thiourea found in reactants, check if it's used for heterocycle formation
            if has_thiourea_in_reactants:
                try:
                    # Check if product has any heterocyclic structure
                    heterocycle_rings = [
                        "benzothiazole",
                        "thiazole",
                        "benzimidazole",
                        "benzoxazole",
                        "oxazole",
                        "imidazole",
                        "triazole",
                        "tetrazole",
                    ]

                    # Check if product contains a heterocycle
                    has_heterocycle = False
                    for ring in heterocycle_rings:
                        if checker.check_ring(ring, product_smiles):
                            print(f"Found heterocycle in product: {ring}")
                            has_heterocycle = True
                            break

                    # Check if thiourea is not present in the product (consumed)
                    thiourea_consumed = not checker.check_fg("Thiourea", product_smiles)

                    # Check if this is a heterocycle formation reaction
                    is_heterocycle_formation = False
                    heterocycle_formation_reactions = [
                        "benzothiazole",
                        "benzimidazole_derivatives_aldehyde",
                        "benzoxazole",
                        "thiazole",
                    ]

                    for rxn_type in heterocycle_formation_reactions:
                        if checker.check_reaction(rxn_type, rsmi):
                            print(f"Found heterocycle formation reaction: {rxn_type}")
                            is_heterocycle_formation = True
                            break

                    # If we have a heterocycle in the product, thiourea is consumed,
                    # and it's a heterocycle formation reaction, then thiourea was used for heterocycle formation
                    if (
                        has_heterocycle
                        and thiourea_consumed
                        and is_heterocycle_formation
                    ):
                        thiourea_used_for_heterocycle = True
                        print("Confirmed thiourea used for heterocycle formation")

                    # If we can't confirm the reaction type but we see thiourea consumed and heterocycle formed,
                    # we'll still consider it a match
                    elif has_heterocycle and thiourea_consumed:
                        thiourea_used_for_heterocycle = True
                        print(
                            "Thiourea likely used for heterocycle formation (reaction type not confirmed)"
                        )

                except Exception as e:
                    print(f"Error checking heterocycle formation: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if thiourea was used as a precursor for heterocycle formation
    return found_thiourea_intermediate and thiourea_used_for_heterocycle
