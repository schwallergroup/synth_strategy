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
    Detects if the synthesis route involves activation of an alcohol to a bromide.
    """
    found_alcohol_to_bromide = False

    def dfs(node, depth=0):
        nonlocal found_alcohol_to_bromide

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["rsmi"]

                # Check for alcohol to bromide conversion reactions
                if (
                    checker.check_reaction("Alkyl bromides from alcohols", rsmi)
                    or checker.check_reaction("PBr3 and alcohol to alkyl bromide", rsmi)
                    or checker.check_reaction("Appel reaction", rsmi)
                ):

                    # Verify that this reaction actually converts an alcohol to a bromide
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    # Check if any reactant has an alcohol group
                    has_alcohol_reactant = any(
                        checker.check_fg("Primary alcohol", r)
                        or checker.check_fg("Secondary alcohol", r)
                        or checker.check_fg("Tertiary alcohol", r)
                        or checker.check_fg("Aromatic alcohol", r)
                        for r in reactants
                    )

                    # Check if product has a primary halide (bromide)
                    has_bromide_product = checker.check_fg("Primary halide", product)

                    if has_alcohol_reactant and has_bromide_product:
                        found_alcohol_to_bromide = True
                        print(
                            f"Found alcohol to bromide activation at depth {depth}: {rsmi}"
                        )
            except Exception as e:
                print(f"Error processing reaction node for alcohol activation: {e}")

        # Recursively process children
        for child in node.get("children", []):
            dfs(child, depth + 1)

    # Start DFS traversal
    dfs(route)

    print(f"Alcohol to bromide activation: {found_alcohol_to_bromide}")
    return found_alcohol_to_bromide
