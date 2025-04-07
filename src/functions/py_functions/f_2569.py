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
    Detects ester hydrolysis as part of the synthetic strategy.
    """
    # Track if we found the pattern
    found_pattern = False

    def dfs_traverse(node, depth=0):
        nonlocal found_pattern

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]

            # Check for all possible hydrolysis reaction types
            hydrolysis_reactions = [
                "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters",
                "Ester saponification (methyl deprotection)",
                "Ester saponification (alkyl deprotection)",
                "COOH ethyl deprotection",
            ]

            if any(checker.check_reaction(rxn, rsmi) for rxn in hydrolysis_reactions):
                found_pattern = True
                print(f"Found ester hydrolysis reaction at depth {depth}: {rsmi}")
                return

            # If specific reaction checks fail, check for functional group changes
            try:
                reactants = rsmi.split(">")[0].split(".")
                products = rsmi.split(">")[2].split(".")

                # Define functional groups to check
                alcohol_types = [
                    "Primary alcohol",
                    "Secondary alcohol",
                    "Tertiary alcohol",
                    "Aromatic alcohol",
                ]
                ester_types = ["Ester", "Carbo-thioester"]

                # Check for functional groups in reactants and products
                has_ester_in_reactants = any(
                    any(checker.check_fg(est, r) for est in ester_types)
                    for r in reactants
                )
                has_carboxylic_acid_in_reactants = any(
                    checker.check_fg("Carboxylic acid", r) for r in reactants
                )
                has_alcohol_in_reactants = any(
                    any(checker.check_fg(alc, r) for alc in alcohol_types)
                    for r in reactants
                )

                has_ester_in_products = any(
                    any(checker.check_fg(est, p) for est in ester_types)
                    for p in products
                )
                has_carboxylic_acid_in_products = any(
                    checker.check_fg("Carboxylic acid", p) for p in products
                )
                has_alcohol_in_products = any(
                    any(checker.check_fg(alc, p) for alc in alcohol_types)
                    for p in products
                )

                # Forward direction: ester → carboxylic acid
                if has_ester_in_reactants and has_carboxylic_acid_in_products:
                    found_pattern = True
                    print(
                        f"Found ester hydrolysis pattern (forward) at depth {depth}: {rsmi}"
                    )
                    return

                # Retrosynthetic direction: carboxylic acid → ester
                # This is the reverse of hydrolysis (esterification) but in retrosynthesis
                # it represents a hydrolysis step in the forward synthesis
                if has_carboxylic_acid_in_reactants and has_ester_in_products:
                    found_pattern = True
                    print(
                        f"Found ester hydrolysis pattern (retrosynthetic) at depth {depth}: {rsmi}"
                    )
                    return

            except Exception as e:
                print(f"Error analyzing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Ester hydrolysis in synthesis: {found_pattern}")
    return found_pattern
