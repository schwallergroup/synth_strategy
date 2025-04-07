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
    This function detects if the synthesis route involves ketone reduction
    in the late stages of the synthesis.

    In retrosynthetic analysis, this means we're looking for reactions where
    a ketone is transformed into a secondary alcohol (forward: alcohol → ketone).
    """
    has_late_ketone_reduction = False
    max_depth_for_late_stage = 1  # Define late stage as depth <= 1

    def dfs_traverse(node, depth=0):
        nonlocal has_late_ketone_reduction

        if node["type"] == "reaction" and depth <= max_depth_for_late_stage:
            try:
                if "rsmi" in node.get("metadata", {}):
                    rsmi = node["metadata"]["rsmi"]
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    print(f"Checking reaction at depth {depth}: {rsmi}")

                    # Method 1: Check directly for reduction reaction types
                    if checker.check_reaction(
                        "Reduction of ketone to secondary alcohol", rsmi
                    ):
                        has_late_ketone_reduction = True
                        print(
                            f"Ketone reduction detected at depth {depth} using reaction checker"
                        )
                        return

                    # Method 2: Check for oxidation reactions (which in retrosynthesis represent reductions)
                    if checker.check_reaction(
                        "Oxidation of alcohols to aldehydes and ketones", rsmi
                    ) or checker.check_reaction(
                        "Oxidation or Dehydrogenation of Alcohols to Aldehydes and Ketones",
                        rsmi,
                    ):
                        # In retrosynthesis, an oxidation reaction is viewed as a reduction
                        # Check if we have secondary alcohol in reactants and ketone in product
                        for reactant in reactants:
                            has_secondary_alcohol = checker.check_fg(
                                "Secondary alcohol", reactant
                            )
                            has_ketone = checker.check_fg("Ketone", product)

                            if has_secondary_alcohol and has_ketone:
                                has_late_ketone_reduction = True
                                print(
                                    f"Ketone reduction detected at depth {depth} (oxidation reaction in forward direction)"
                                )
                                return

                    # Method 3: General check for functional group transformation
                    # In retrosynthesis, we're looking for ketone in product and secondary alcohol in reactants
                    has_ketone_in_product = checker.check_fg("Ketone", product)

                    if has_ketone_in_product:
                        for reactant in reactants:
                            if checker.check_fg("Secondary alcohol", reactant):
                                # This pattern suggests alcohol oxidation (reduction in retrosynthesis)
                                has_late_ketone_reduction = True
                                print(
                                    f"Ketone reduction detected at depth {depth} using FG transformation check"
                                )
                                return
            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)

    if not has_late_ketone_reduction:
        print("No late-stage ketone reduction detected in the synthesis route")

    return has_late_ketone_reduction
