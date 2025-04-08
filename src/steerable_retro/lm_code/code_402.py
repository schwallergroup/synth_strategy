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
    This function detects if the route employs a strategy involving both oxidation and reduction
    of sulfur-containing functional groups (thioether ↔ sulfone cycling).
    """
    # Track sulfur oxidation and reduction events
    oxidation_events = []  # (depth, reactant_smiles, product_smiles)
    reduction_events = []  # (depth, reactant_smiles, product_smiles)

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            # Extract reactants and product
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for sulfur oxidation and reduction based on functional groups
                has_oxidation = False
                has_reduction = False

                # Check for oxidation patterns
                if (
                    any(checker.check_fg("Monosulfide", r) for r in reactants)
                    and checker.check_fg("Sulfoxide", product)
                    and not any(checker.check_fg("Sulfoxide", r) for r in reactants)
                ):
                    print(f"Found sulfur oxidation (thioether → sulfoxide) at depth {depth}")
                    oxidation_events.append((depth, reactants, product))
                    has_oxidation = True

                if (
                    any(checker.check_fg("Monosulfide", r) for r in reactants)
                    and checker.check_fg("Sulfone", product)
                    and not any(checker.check_fg("Sulfone", r) for r in reactants)
                ):
                    print(f"Found sulfur oxidation (thioether → sulfone) at depth {depth}")
                    oxidation_events.append((depth, reactants, product))
                    has_oxidation = True

                if (
                    any(checker.check_fg("Sulfoxide", r) for r in reactants)
                    and checker.check_fg("Sulfone", product)
                    and not any(checker.check_fg("Sulfone", r) for r in reactants)
                ):
                    print(f"Found sulfur oxidation (sulfoxide → sulfone) at depth {depth}")
                    oxidation_events.append((depth, reactants, product))
                    has_oxidation = True

                # Check for reduction patterns
                if (
                    any(checker.check_fg("Sulfone", r) for r in reactants)
                    and checker.check_fg("Sulfoxide", product)
                    and not checker.check_fg("Sulfone", product)
                ):
                    print(f"Found sulfur reduction (sulfone → sulfoxide) at depth {depth}")
                    reduction_events.append((depth, reactants, product))
                    has_reduction = True

                if (
                    any(checker.check_fg("Sulfone", r) for r in reactants)
                    and checker.check_fg("Monosulfide", product)
                    and not checker.check_fg("Sulfone", product)
                ):
                    print(f"Found sulfur reduction (sulfone → thioether) at depth {depth}")
                    reduction_events.append((depth, reactants, product))
                    has_reduction = True

                if (
                    any(checker.check_fg("Sulfoxide", r) for r in reactants)
                    and checker.check_fg("Monosulfide", product)
                    and not checker.check_fg("Sulfoxide", product)
                ):
                    print(f"Found sulfur reduction (sulfoxide → thioether) at depth {depth}")
                    reduction_events.append((depth, reactants, product))
                    has_reduction = True

                # Check for specific reaction types
                rxn_smiles = node["metadata"].get("smiles", "")
                if rxn_smiles:
                    # Check for oxidation reactions
                    if (
                        checker.check_reaction("Sulfanyl to sulfinyl", rxn_smiles)
                        or checker.check_reaction("Sulfanyl to sulfinyl_peroxide", rxn_smiles)
                        or checker.check_reaction("Sulfanyl to sulfinyl_H2O2", rxn_smiles)
                        or checker.check_reaction("Sulfanyl to sulfinyl_H2O", rxn_smiles)
                        or checker.check_reaction("Sulfanyl to sulfinyl_SO3-", rxn_smiles)
                        or checker.check_reaction("Sulfanyl to sulfinyl_sulfonyl", rxn_smiles)
                        or checker.check_reaction("Sulfanyl to sulfinyl_MeOH", rxn_smiles)
                        or checker.check_reaction("Sulfanyl to sulfinyl_COO", rxn_smiles)
                    ):
                        print(f"Found sulfur oxidation reaction at depth {depth}")
                        if not has_oxidation:  # Avoid duplicates
                            oxidation_events.append((depth, reactants, product))

                    # Check for reduction reactions - look for thia-Michael addition in reverse
                    # (since we're traversing retrosynthetically)
                    if checker.check_reaction("thia-Michael addition", rxn_smiles):
                        # In retrosynthesis, thia-Michael addition can represent a reduction
                        # when traversing from product to reactants
                        print(f"Found potential sulfur reduction reaction at depth {depth}")
                        if not has_reduction:  # Avoid duplicates
                            reduction_events.append((depth, reactants, product))

                    # Check for cleavage of sulfones and sulfoxides (reduction)
                    if checker.check_reaction("Cleavage of sulfons and sulfoxides", rxn_smiles):
                        print(f"Found sulfur reduction reaction at depth {depth}")
                        if not has_reduction:  # Avoid duplicates
                            reduction_events.append((depth, reactants, product))

            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if both oxidation and reduction events were found
    has_oxidation = len(oxidation_events) > 0
    has_reduction = len(reduction_events) > 0

    # For a true redox cycling strategy, we need both oxidation and reduction
    # of sulfur-containing functional groups
    has_redox_cycling = has_oxidation and has_reduction

    print(f"Oxidation events: {len(oxidation_events)}, Reduction events: {len(reduction_events)}")
    print(f"Has redox cycling: {has_redox_cycling}")

    # If we have oxidation but no reduction, check if any of the reactions might be
    # misclassified reductions (some reactions can be hard to classify)
    if has_oxidation and not has_reduction:
        # Look for any reaction that might involve sulfur and could potentially be a reduction
        def check_potential_reduction(node):
            if node["type"] == "reaction":
                try:
                    rsmi = node["metadata"]["rsmi"]
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    # Check if product contains sulfur and reactants contain oxidizing agents
                    if (
                        (any("S" in r for r in reactants) and "S" in product)
                        or (
                            any(checker.check_fg("Monosulfide", r) for r in reactants)
                            and checker.check_fg("Monosulfide", product)
                        )
                        or (
                            any(checker.check_fg("Sulfoxide", r) for r in reactants)
                            and checker.check_fg("Monosulfide", product)
                        )
                        or (
                            any(checker.check_fg("Sulfone", r) for r in reactants)
                            and checker.check_fg("Monosulfide", product)
                        )
                        or (
                            any(checker.check_fg("Sulfone", r) for r in reactants)
                            and checker.check_fg("Sulfoxide", product)
                        )
                    ):
                        print(f"Found potential unlabeled sulfur reduction")
                        return True
                except Exception:
                    pass

            for child in node.get("children", []):
                if check_potential_reduction(child):
                    return True
            return False

        # If we find a potential reduction, consider it a redox cycling strategy
        if check_potential_reduction(route):
            print("Detected potential unlabeled sulfur reduction - considering as redox cycling")
            has_redox_cycling = True

    return has_redox_cycling
