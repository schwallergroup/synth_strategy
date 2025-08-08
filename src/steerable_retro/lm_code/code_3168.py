#!/bin/python

"""LM-defined function for strategy description."""

from rdkit.Chem import AllChem, rdFMCS
import copy
from collections import deque
import rdkit.Chem as Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdChemReactions
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
import rdkit.Chem.rdFMCS
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem import AllChem, Descriptors, Lipinski
from rdkit.Chem import rdmolops
import re
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import AllChem, Descriptors
import traceback
import rdkit
from collections import Counter
from steerable_retro.utils.check import Check
from steerable_retro.utils import fuzzy_dict, check

root_data = "/home/dparm/steerable_retro/data"

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
    Detects a linear synthesis strategy that includes a protection step and maintains
    the protecting group throughout the synthesis.
    """
    # Track protection and deprotection steps
    protection_step_found = False
    deprotection_step_found = False

    # Track reactions with protected intermediates
    reactions_with_protection = []
    all_reactions = []

    # Track depth of protection and deprotection
    protection_depth = -1
    deprotection_depth = -1

    def dfs_traverse(node, depth):
        nonlocal protection_step_found, deprotection_step_found
        nonlocal protection_depth, deprotection_depth

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0]
            product = rsmi.split(">")[-1]

            all_reactions.append((depth, rsmi))

            # Check for protection step (ketone/aldehyde to acetal/ketal)
            if (
                (checker.check_fg("Ketone", reactants) or checker.check_fg("Aldehyde", reactants))
                and checker.check_ring("dioxolane", product)
                and checker.check_reaction("Aldehyde or ketone acetalization", rsmi)
            ) or (
                checker.check_fg("Diol", reactants)
                and checker.check_ring("dioxolane", product)
                and checker.check_reaction("Diol acetalization", rsmi)
            ):
                protection_step_found = True
                protection_depth = depth
                reactions_with_protection.append((depth, rsmi))
                print(f"Found protection step at depth {depth}: {rsmi}")

            # Check for deprotection step (acetal/ketal to ketone/aldehyde)
            elif (
                checker.check_ring("dioxolane", reactants)
                and (checker.check_fg("Ketone", product) or checker.check_fg("Aldehyde", product))
                and (
                    checker.check_reaction("Acetal hydrolysis to aldehyde", rsmi)
                    or checker.check_reaction("Ketal hydrolysis to ketone", rsmi)
                )
            ):
                deprotection_step_found = True
                deprotection_depth = depth
                print(f"Found deprotection step at depth {depth}: {rsmi}")

            # Check if dioxolane is maintained in intermediate steps
            elif checker.check_ring("dioxolane", reactants) and checker.check_ring(
                "dioxolane", product
            ):
                reactions_with_protection.append((depth, rsmi))
                # If we haven't found an explicit protection step yet, consider this as evidence
                if not protection_step_found:
                    protection_step_found = True
                    protection_depth = depth
                print(f"Dioxolane protecting group maintained at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route, 0)

    # Sort reactions by depth to analyze sequence
    all_reactions.sort(key=lambda x: x[0])
    reactions_with_protection.sort(key=lambda x: x[0])

    # Strategy criteria:
    # 1. Protection is maintained across multiple reactions
    # 2. If deprotection occurs, it should be at a later stage than protection

    has_protection_strategy = False

    # If we have at least two reactions with protected intermediates, that's a protection strategy
    if len(reactions_with_protection) >= 2:
        has_protection_strategy = True
    # Or if we have explicit protection and deprotection in the right order
    elif protection_step_found and deprotection_step_found:
        has_protection_strategy = protection_depth < deprotection_depth

    print(f"Linear synthesis with protection strategy detected: {has_protection_strategy}")
    print(f"Protection depth: {protection_depth}, Deprotection depth: {deprotection_depth}")
    print(
        f"Total reactions: {len(all_reactions)}, Protected reactions: {len(reactions_with_protection)}"
    )

    return has_protection_strategy
