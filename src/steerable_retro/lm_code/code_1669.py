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
    Detects a strategy where nitrogen-containing groups are incorporated in the late stages
    of synthesis through nucleophilic aromatic substitution of halogens.
    """
    # Track key features
    n_incorporations = 0
    n_incorporation_depths = []

    # Define nitrogen-containing functional groups to check
    n_functional_groups = [
        "Primary amine",
        "Secondary amine",
        "Tertiary amine",
        "Aniline",
        "Azide",
        "Nitro group",
        "Nitrile",
        "Primary amide",
        "Secondary amide",
        "Tertiary amide",
        "Hydrazine",
        "Hydrazone",
    ]

    # Define N-arylation reaction types
    n_arylation_reactions = [
        "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
        "Goldberg coupling",
        "Goldberg coupling aryl amine-aryl chloride",
        "Goldberg coupling aryl amide-aryl chloride",
        "Ullmann-Goldberg Substitution amine",
        "Buchwald-Hartwig",
        "Ullmann condensation",
        "N-arylation_heterocycles",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal n_incorporations, n_incorporation_depths

        if node["type"] == "reaction":
            # Extract reactants and product
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check for nitrogen incorporation via halogen displacement
                n_nucleophile = False
                n_nucleophile_smiles = None
                halogenated_aromatic = False
                halogenated_reactant = None

                # Check for nitrogen nucleophiles in reactants
                for reactant_smiles in reactants_smiles:
                    for n_fg in n_functional_groups:
                        if checker.check_fg(n_fg, reactant_smiles):
                            n_nucleophile = True
                            n_nucleophile_smiles = reactant_smiles
                            print(f"Found nitrogen nucleophile ({n_fg}): {reactant_smiles}")
                            break

                    # Check for halogenated aromatic
                    if checker.check_fg("Aromatic halide", reactant_smiles):
                        halogenated_aromatic = True
                        halogenated_reactant = reactant_smiles
                        print(f"Found aromatic halide: {reactant_smiles}")

                # Check if this is an N-arylation reaction
                is_n_arylation = False
                for rxn_type in n_arylation_reactions:
                    if checker.check_reaction(rxn_type, rsmi):
                        is_n_arylation = True
                        print(f"Detected {rxn_type} reaction")
                        break

                # If no specific reaction type was found, check for general nucleophilic aromatic substitution
                if not is_n_arylation and n_nucleophile and halogenated_aromatic:
                    # This is a fallback check for nucleophilic aromatic substitution
                    print("Checking for general nucleophilic aromatic substitution")

                    # Check if the product has a new C-N bond where the halogen was
                    if any(checker.check_fg(n_fg, product_smiles) for n_fg in n_functional_groups):
                        # If the halogenated reactant didn't have nitrogen but the product does
                        if not any(
                            checker.check_fg(n_fg, halogenated_reactant)
                            for n_fg in n_functional_groups
                        ):
                            is_n_arylation = True
                            print("Detected general nucleophilic aromatic substitution")

                # If conditions are met and depth is low (late stage)
                if n_nucleophile and halogenated_aromatic and is_n_arylation and depth <= 3:
                    print(f"Potential nitrogen incorporation at depth {depth}")

                    # Check if the halogenated reactant doesn't already have the nitrogen group
                    if not any(
                        checker.check_fg(n_fg, halogenated_reactant) for n_fg in n_functional_groups
                    ):
                        # Verify nitrogen is incorporated in product
                        if any(
                            checker.check_fg(n_fg, product_smiles) for n_fg in n_functional_groups
                        ):
                            # Additional check: ensure the nitrogen nucleophile and halogenated aromatic
                            # are different reactants (not the same molecule)
                            if n_nucleophile_smiles != halogenated_reactant:
                                n_incorporations += 1
                                n_incorporation_depths.append(depth)
                                print(f"Confirmed nitrogen incorporation at depth {depth}")
                                print(f"Halogenated reactant: {halogenated_reactant}")
                                print(f"Nitrogen nucleophile: {n_nucleophile_smiles}")
                                print(f"Product with nitrogen: {product_smiles}")
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if the strategy is present
    # We need at least 1 nitrogen incorporation in the late stages (depth <= 3)
    strategy_present = n_incorporations >= 1

    print(f"Nitrogen incorporations: {n_incorporations}")
    print(f"Incorporation depths: {n_incorporation_depths}")
    print(f"Strategy present: {strategy_present}")

    return strategy_present
