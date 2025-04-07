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
    Detects a synthesis where nitro reduction is used as the final step or penultimate step.
    """
    has_nitro_reduction = False
    nitro_reduction_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal has_nitro_reduction, nitro_reduction_depth

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for nitro reduction using the checker function
                if checker.check_reaction("Reduction of nitro groups to amines", rsmi):
                    has_nitro_reduction = True
                    nitro_reduction_depth = depth
                    print(f"Found nitro reduction at depth {depth}, rsmi: {rsmi}")
                # Fallback check in case the reaction checker doesn't catch it
                elif any(checker.check_fg("Nitro group", reactant) for reactant in reactants) and (
                    checker.check_fg("Primary amine", product)
                    or checker.check_fg("Aniline", product)
                ):
                    # Verify it's actually a reduction by checking for nitro group disappearance
                    if not checker.check_fg("Nitro group", product):
                        has_nitro_reduction = True
                        nitro_reduction_depth = depth
                        print(
                            f"Found nitro reduction (fallback method) at depth {depth}, rsmi: {rsmi}"
                        )
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Consider nitro reduction as late-stage if it occurs at depth 0 or 1
    # Depth 0 is the final step, depth 1 is the penultimate step
    strategy_detected = has_nitro_reduction and nitro_reduction_depth <= 1

    print(f"Late-stage nitro reduction strategy detected: {strategy_detected}")
    print(f"Nitro reduction found: {has_nitro_reduction}, at depth: {nitro_reduction_depth}")
    return strategy_detected
