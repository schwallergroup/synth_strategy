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
    Detects nitrile hydrolysis as the final step in the synthesis.
    """
    found_nitrile_hydrolysis = False

    def dfs_traverse(node, depth=0):
        nonlocal found_nitrile_hydrolysis

        if node["type"] == "reaction" and depth <= 1:  # Check final or penultimate step
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if this is a nitrile hydrolysis reaction
                if checker.check_reaction(
                    "Oxidation of nitrile to carboxylic acid", rsmi
                ):
                    print(f"Found nitrile hydrolysis reaction at depth {depth}: {rsmi}")
                    found_nitrile_hydrolysis = True
                else:
                    # Double-check by looking at functional groups
                    reactant_has_nitrile = any(
                        checker.check_fg("Nitrile", r) for r in reactants
                    )
                    product_has_carboxylic_acid = checker.check_fg(
                        "Carboxylic acid", product
                    )
                    nitrile_in_product = checker.check_fg("Nitrile", product)

                    if (
                        reactant_has_nitrile
                        and product_has_carboxylic_acid
                        and not nitrile_in_product
                    ):
                        print(
                            f"Found nitrile hydrolysis by functional group check at depth {depth}: {rsmi}"
                        )
                        found_nitrile_hydrolysis = True
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Late-stage nitrile hydrolysis detected: {found_nitrile_hydrolysis}")
    return found_nitrile_hydrolysis
