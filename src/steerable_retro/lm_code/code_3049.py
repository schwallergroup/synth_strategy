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
    This function detects late-stage ketone reduction (C=O to C-OH).
    Late stage is defined as occurring at depth 0 or 1.

    In retrosynthetic analysis, we look for the reverse transformation:
    secondary alcohol in reactants being converted to ketone in product.
    """
    late_stage_reduction_found = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_reduction_found

        if node["type"] == "reaction" and depth <= 1:
            # Extract reactants and product
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                print(f"Examining reaction at depth {depth}: {rsmi}")
                print(f"  Reactants: {reactants_smiles}")
                print(f"  Product: {product_smiles}")

                # Check for alcohol oxidation in retrosynthesis
                oxidation_reaction = checker.check_reaction(
                    "Oxidation or Dehydrogenation of Alcohols to Aldehydes and Ketones", rsmi
                )

                if not oxidation_reaction:
                    # Try alternative reaction types that might involve alcohol oxidation
                    oxidation_reaction = checker.check_reaction(
                        "Oxidation of secondary alcohol to ketone", rsmi
                    )

                print(f"  Is oxidation reaction: {oxidation_reaction}")

                if oxidation_reaction:
                    # In retrosynthesis: alcohol in reactants, ketone in product
                    has_secondary_alcohol = any(
                        checker.check_fg("Secondary alcohol", reactant)
                        for reactant in reactants_smiles
                    )
                    has_ketone = checker.check_fg("Ketone", product_smiles)

                    print(f"  Ketone in product: {has_ketone}")
                    print(f"  Secondary alcohol in reactants: {has_secondary_alcohol}")

                    if has_ketone and has_secondary_alcohol:
                        print(f"Found late-stage ketone reduction at depth {depth}")
                        print(f"Reaction SMILES: {rsmi}")
                        late_stage_reduction_found = True

                # If we didn't find the pattern through reaction type checking,
                # try direct functional group analysis
                if not late_stage_reduction_found:
                    has_secondary_alcohol = any(
                        checker.check_fg("Secondary alcohol", reactant)
                        for reactant in reactants_smiles
                    )
                    has_ketone = checker.check_fg("Ketone", product_smiles)

                    # Check if this is a simple alcohol to ketone transformation
                    if has_secondary_alcohol and has_ketone and depth <= 1:
                        # Additional check to ensure it's not a different reaction type
                        if not any(
                            checker.check_reaction(rxn, rsmi)
                            for rxn in [
                                "Acetal hydrolysis to ketone",
                                "Ketal hydrolysis to ketone",
                                "Grignard from ketone to alcohol",
                            ]
                        ):
                            print(
                                f"Found late-stage ketone reduction through FG analysis at depth {depth}"
                            )
                            print(f"Reaction SMILES: {rsmi}")
                            late_stage_reduction_found = True

            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return late_stage_reduction_found
