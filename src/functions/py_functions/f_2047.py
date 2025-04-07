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
    Detects if the synthesis route involves an early-stage conversion of
    a primary amine to an azide group (typically depth > 1).

    Note: In retrosynthesis, we're looking for azide → amine transformations
    in the reaction SMILES (which are in forward direction).
    """
    found_amine_to_azide = False

    def dfs_traverse(node, depth=0):
        nonlocal found_amine_to_azide

        if node["type"] == "reaction" and depth > 1:  # Early stage (depth > 1)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Checking reaction at depth {depth}: {rsmi}")

                # In retrosynthesis, we're looking at forward reactions
                # where azide is in reactants and amine is in product

                # Check if this is an "Amine to azide" reaction
                if checker.check_reaction("Amine to azide", rsmi):
                    print(f"Found 'Amine to azide' reaction at depth {depth}")
                    found_amine_to_azide = True
                else:
                    # Alternative check: verify azide in reactants and primary amine in product
                    has_azide = any(
                        checker.check_fg("Azide", r) for r in reactants if r
                    )
                    has_primary_amine = checker.check_fg("Primary amine", product)

                    print(
                        f"Azide in reactants: {has_azide}, Primary amine in product: {has_primary_amine}"
                    )

                    if has_azide and has_primary_amine:
                        print(
                            f"Found azide in reactants and primary amine in product at depth {depth}"
                        )
                        found_amine_to_azide = True

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    return found_amine_to_azide
