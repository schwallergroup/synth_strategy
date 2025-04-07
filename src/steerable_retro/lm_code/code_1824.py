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
    This function detects a synthetic strategy involving oxime intermediates
    leading to isoxazole formation, followed by late-stage amide coupling.
    """
    # Track the sequence of reactions and intermediates
    reaction_sequence = []

    # Track if we've found isoxazole in the synthesis
    isoxazole_found = False

    def dfs_traverse(node, depth=0):
        nonlocal isoxazole_found, reaction_sequence

        # Check for isoxazole in molecule nodes
        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            if checker.check_ring("isoxazole", mol_smiles):
                isoxazole_found = True
                print(f"Found isoxazole ring at depth {depth}: {mol_smiles}")

        # Check reactions
        elif node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check for oxime formation
                if checker.check_fg("Oxime", product_smiles) and any(
                    checker.check_fg("Aldehyde", r) for r in reactants_smiles
                ):
                    reaction_sequence.append(("oxime_formation", depth))
                    print(f"Detected oxime formation at depth {depth}")

                # Check for chloro-oxime formation
                if (
                    any(checker.check_fg("Oxime", r) for r in reactants_smiles)
                    and checker.check_fg("Substituted imine", product_smiles)
                    and any("Cl" in r for r in reactants_smiles)
                ):
                    reaction_sequence.append(("chloro_oxime", depth))
                    print(f"Detected chloro-oxime formation at depth {depth}")

                # Check for isoxazole formation
                if checker.check_ring("isoxazole", product_smiles) and not any(
                    checker.check_ring("isoxazole", r) for r in reactants_smiles
                ):
                    reaction_sequence.append(("isoxazole_formation", depth))
                    print(f"Detected isoxazole formation at depth {depth}")

                # Check for ester hydrolysis
                if checker.check_reaction(
                    "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters", rsmi
                ):
                    reaction_sequence.append(("ester_hydrolysis", depth))
                    print(f"Detected ester hydrolysis at depth {depth}")

                # Check for amide coupling
                if any(checker.check_fg("Carboxylic acid", r) for r in reactants_smiles) and (
                    checker.check_fg("Primary amide", product_smiles)
                    or checker.check_fg("Secondary amide", product_smiles)
                    or checker.check_fg("Tertiary amide", product_smiles)
                ):
                    reaction_sequence.append(("amide_coupling", depth))
                    print(f"Detected amide coupling at depth {depth}")

                # Alternative check for amide coupling using reaction type
                if (
                    checker.check_reaction(
                        "Acylation of Nitrogen Nucleophiles by Carboxylic Acids", rsmi
                    )
                    or checker.check_reaction("Carboxylic acid with primary amine to amide", rsmi)
                    or checker.check_reaction(
                        "Acyl chloride with primary amine to amide (Schotten-Baumann)", rsmi
                    )
                ):
                    if not any(rxn == "amide_coupling" for rxn, _ in reaction_sequence):
                        reaction_sequence.append(("amide_coupling", depth))
                        print(f"Detected amide coupling reaction at depth {depth}")

            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Sort reactions by depth to get the sequence
    reaction_sequence.sort(key=lambda x: x[1], reverse=True)
    reaction_types = [r[0] for r in reaction_sequence]

    print(f"Detected reaction sequence: {reaction_types}")

    # Check if we have the complete strategy in the correct order
    # The expected sequence is: oxime_formation -> chloro_oxime -> isoxazole_formation -> ester_hydrolysis -> amide_coupling
    expected_sequence = [
        "oxime_formation",
        "chloro_oxime",
        "isoxazole_formation",
        "ester_hydrolysis",
        "amide_coupling",
    ]

    # Check if all expected reactions are present
    all_reactions_present = all(rxn in reaction_types for rxn in expected_sequence)

    # Check if the reactions are in the correct order
    correct_order = True
    for i in range(len(expected_sequence) - 1):
        if expected_sequence[i] in reaction_types and expected_sequence[i + 1] in reaction_types:
            idx1 = reaction_types.index(expected_sequence[i])
            idx2 = reaction_types.index(expected_sequence[i + 1])
            if idx1 > idx2:  # Earlier reaction should have higher depth
                correct_order = False
                break

    # If isoxazole is found but not all reactions are detected, we can still consider it a partial match
    # This handles cases where some intermediate steps might be missing in the route
    if (
        isoxazole_found
        and "isoxazole_formation" in reaction_types
        and "amide_coupling" in reaction_types
    ):
        partial_match = True
        print("Detected partial isoxazole via oxime strategy (key steps present)")
    else:
        partial_match = False

    strategy_present = all_reactions_present and correct_order and isoxazole_found

    if strategy_present:
        print("Detected complete isoxazole via oxime strategy")
    else:
        print("Incomplete strategy detection:")
        print(f"  All reactions present: {all_reactions_present}")
        print(f"  Correct reaction order: {correct_order}")
        print(f"  Isoxazole found: {isoxazole_found}")
        print(f"  Reaction sequence: {reaction_types}")

    # Return true if we have either a complete match or a partial match with key steps
    return strategy_present or partial_match
