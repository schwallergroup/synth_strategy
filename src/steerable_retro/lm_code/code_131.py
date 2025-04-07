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
    Detects if the synthesis route includes late-stage esterification via acid chloride activation.
    """
    # Track if we found acid chloride formation and subsequent esterification
    found_acid_chloride = False
    found_esterification = False
    acid_chloride_depth = float("inf")
    max_depth = get_max_depth(route)

    def dfs(node, depth=0):
        nonlocal found_acid_chloride, found_esterification, acid_chloride_depth

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]

            # Check for acid chloride formation
            if checker.check_reaction("Acyl chlorides from alcohols", rsmi):
                found_acid_chloride = True
                acid_chloride_depth = depth
                print(f"Found acid chloride formation at depth {depth}")

            # Also check for acyl halide in product
            product = rsmi.split(">")[-1]
            if checker.check_fg("Acyl halide", product):
                found_acid_chloride = True
                acid_chloride_depth = min(acid_chloride_depth, depth)
                print(f"Found acyl halide formation at depth {depth}")

            # Check for esterification using acid chloride
            if checker.check_reaction("Schotten-Baumann to ester", rsmi) or checker.check_reaction(
                "Esterification of Carboxylic Acids", rsmi
            ):
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if one reactant has acid chloride and product has ester
                has_acid_chloride = any(checker.check_fg("Acyl halide", r) for r in reactants)
                has_ester = checker.check_fg("Ester", product)

                if has_acid_chloride and has_ester:
                    found_esterification = True
                    print(f"Found esterification via acid chloride at depth {depth}")

        # Continue DFS traversal
        for child in node.get("children", []):
            dfs(child, depth + 1)

    # Start DFS traversal
    dfs(route)

    # Ensure esterification is late-stage (low depth relative to max depth)
    # Late-stage is defined as occurring in the first half of the synthesis depth
    late_stage_threshold = max(1, max_depth // 2)
    print(
        f"Acid chloride depth: {acid_chloride_depth}, Late stage threshold: {late_stage_threshold}"
    )

    return (
        found_acid_chloride and found_esterification and acid_chloride_depth <= late_stage_threshold
    )
