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
    This function detects ester reduction to alcohol as part of the synthetic strategy.
    """
    ester_reduction = False

    def dfs_traverse(node):
        nonlocal ester_reduction

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this is an ester reduction reaction directly
            if checker.check_reaction(
                "Reduction of ester to primary alcohol", rsmi
            ) or checker.check_reaction(
                "Reduction of carboxylic acid to primary alcohol", rsmi
            ):
                print(f"Found ester/acid reduction reaction: {rsmi}")
                ester_reduction = True
            else:
                # In retrosynthesis, primary alcohol is in reactants, ester is in product
                has_primary_alcohol = any(
                    checker.check_fg("Primary alcohol", r) for r in reactants
                )
                has_ester = checker.check_fg("Ester", product)

                if has_primary_alcohol and has_ester:
                    print(
                        f"Found potential ester reduction (retrosynthetic): primary alcohol in reactants, ester in product: {rsmi}"
                    )
                    ester_reduction = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(f"Ester reduction strategy detected: {ester_reduction}")

    return ester_reduction
