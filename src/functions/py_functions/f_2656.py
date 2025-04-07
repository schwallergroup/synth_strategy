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
    Checks if the synthesis route contains a nitro reduction as one of the final steps
    """

    def dfs_traverse(node, depth, nitro_reductions):
        if node["type"] == "reaction":
            # Check if this is a nitro reduction reaction
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if any reactant has a nitro group
                has_nitro_reactant = any(
                    checker.check_fg("Nitro group", reactant) for reactant in reactants
                )

                # Check if product has an amine group
                has_amine_product = (
                    checker.check_fg("Primary amine", product)
                    or checker.check_fg("Secondary amine", product)
                    or checker.check_fg("Tertiary amine", product)
                )

                # Check if this is a nitro reduction reaction
                is_nitro_reduction = checker.check_reaction(
                    "Reduction of nitro groups to amines", rsmi
                )

                if (has_nitro_reactant and has_amine_product) or is_nitro_reduction:
                    nitro_reductions.append((depth, rsmi))
            except Exception as e:
                print(f"Error analyzing reaction: {e}")

        # Recursively traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1, nitro_reductions)

    # Start traversal
    nitro_reductions = []
    dfs_traverse(route, 0, nitro_reductions)

    # Sort by depth to find the latest (smallest depth) nitro reduction
    if nitro_reductions:
        nitro_reductions.sort(key=lambda x: x[0])
        latest_nitro_reduction = nitro_reductions[0]

        # Check if it's one of the final steps (depth ≤ 2)
        if latest_nitro_reduction[0] <= 2:
            print(f"Found nitro reduction as final step: {latest_nitro_reduction}")
            return True

    return False
