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
    This function detects a synthetic strategy involving isothiocyanate intermediate,
    thiourea formation, and heterocycle construction with late-stage ring formation.
    """
    # Initialize tracking variables
    isothiocyanate_reactions = []
    thiourea_reactions = []
    heterocycle_formations = []

    # Track the sequence of transformations
    reaction_sequence = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            # Extract reactants and product
            try:
                rsmi = node.get("metadata", {}).get("rsmi", "")
                if not rsmi:
                    return

                parts = rsmi.split(">")
                if len(parts) < 3:
                    return

                reactants_smiles = parts[0].split(".")
                product_smiles = parts[2]

                # Store reaction with depth for sequence analysis
                reaction_data = {
                    "rsmi": rsmi,
                    "depth": depth,
                    "reactants": reactants_smiles,
                    "product": product_smiles,
                }
                reaction_sequence.append(reaction_data)

                # Check for isothiocyanate in reactants
                for reactant in reactants_smiles:
                    if checker.check_fg("Isothiocyanate", reactant):
                        print(f"Found isothiocyanate in reactant: {reactant}")
                        isothiocyanate_reactions.append(reaction_data)

                # Check for thiourea in product
                if checker.check_fg("Thiourea", product_smiles):
                    print(f"Found thiourea in product: {product_smiles}")
                    thiourea_reactions.append(reaction_data)

                # Check for heterocycle formation
                heterocycle_formed = False
                heterocycle_types = [
                    "benzothiazole",
                    "thiazole",
                    "benzimidazole",
                    "imidazole",
                    "oxazole",
                    "benzoxazole",
                    "triazole",
                    "tetrazole",
                ]

                # Check if product contains a heterocycle that reactants don't have
                for ring_type in heterocycle_types:
                    if checker.check_ring(ring_type, product_smiles):
                        # Check if any reactant already has this heterocycle
                        reactant_has_ring = any(
                            checker.check_ring(ring_type, r) for r in reactants_smiles
                        )
                        if not reactant_has_ring:
                            print(
                                f"Found heterocycle formation ({ring_type}) at depth {depth}"
                            )
                            heterocycle_formed = True
                            heterocycle_formations.append(reaction_data)
                            break

            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Process children with increased depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Sort reactions by depth to analyze the sequence
    reaction_sequence.sort(key=lambda x: x["depth"])

    # Check if we have the complete strategy
    has_isothiocyanate = len(isothiocyanate_reactions) > 0
    has_thiourea = len(thiourea_reactions) > 0
    has_heterocycle_formation = len(heterocycle_formations) > 0

    # Check if the sequence is correct (isothiocyanate -> thiourea -> heterocycle)
    correct_sequence = False
    if has_isothiocyanate and has_thiourea and has_heterocycle_formation:
        # Get the minimum depths for each step
        isothiocyanate_depth = min(r["depth"] for r in isothiocyanate_reactions)
        thiourea_depth = min(r["depth"] for r in thiourea_reactions)
        heterocycle_depth = min(r["depth"] for r in heterocycle_formations)

        # Check if the sequence is in the correct order
        # In retrosynthetic direction: heterocycle (lowest depth) -> thiourea -> isothiocyanate (highest depth)
        if heterocycle_depth <= thiourea_depth <= isothiocyanate_depth:
            print(
                f"Correct sequence detected: heterocycle ({heterocycle_depth}) -> thiourea ({thiourea_depth}) -> isothiocyanate ({isothiocyanate_depth})"
            )
            correct_sequence = True

    # Strategy is present if we have all components in the correct sequence
    strategy_present = (
        has_isothiocyanate
        and has_thiourea
        and has_heterocycle_formation
        and correct_sequence
    )

    print(f"Has isothiocyanate: {has_isothiocyanate}")
    print(f"Has thiourea: {has_thiourea}")
    print(f"Has heterocycle formation: {has_heterocycle_formation}")
    print(f"Correct sequence: {correct_sequence}")
    print(f"Isothiocyanate-thiourea heterocycle strategy detected: {strategy_present}")

    return strategy_present
