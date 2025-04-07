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
    Detects if the synthesis includes a nitro reduction to amine sequence.
    """
    found_nitro_reduction = False

    def dfs_traverse(node, depth=0):
        nonlocal found_nitro_reduction

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]

                # Check if this is a nitro reduction reaction using the checker
                if checker.check_reaction("Reduction of nitro groups to amines", rsmi):
                    print(f"Found nitro reduction reaction at depth {depth}: {rsmi}")
                    found_nitro_reduction = True
                else:
                    # Fallback method: check reactants and products manually
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    # Check if product contains primary amine
                    has_amine_in_product = checker.check_fg(
                        "Primary amine", product
                    ) or checker.check_fg("Aniline", product)

                    if has_amine_in_product:
                        # Check if any reactant contains nitro group
                        for reactant in reactants:
                            if checker.check_fg("Nitro group", reactant):
                                print(
                                    f"Found nitro reduction sequence at depth {depth} (manual check)"
                                )
                                print(f"Reaction: {rsmi}")
                                print(f"Reactant with nitro group: {reactant}")
                                print(f"Product with amine: {product}")
                                found_nitro_reduction = True
                                break

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    return found_nitro_reduction
