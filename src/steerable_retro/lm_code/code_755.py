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
    Detects if the synthesis uses aldehyde as a key early intermediate.
    """
    aldehyde_intermediate_detected = False
    min_depth_for_early_stage = 4  # Define what "early stage" means

    def dfs_traverse(node, depth=0):
        nonlocal aldehyde_intermediate_detected

        if node["type"] == "reaction" and depth >= min_depth_for_early_stage:
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0]
                product_smiles = rsmi.split(">")[-1]

                # Check if product contains aldehyde
                product_has_aldehyde = checker.check_fg("Aldehyde", product_smiles)

                # Check if any reactant contains aldehyde
                reactants = reactants_smiles.split(".")
                reactant_has_aldehyde = any(checker.check_fg("Aldehyde", r) for r in reactants)

                # Aldehyde is an intermediate if it's either formed or consumed in the reaction
                if (product_has_aldehyde and not reactant_has_aldehyde) or (
                    not product_has_aldehyde and reactant_has_aldehyde
                ):
                    print(f"Aldehyde intermediate detected at depth {depth}")
                    print(f"Reaction SMILES: {rsmi}")
                    if product_has_aldehyde and not reactant_has_aldehyde:
                        print("Aldehyde is formed in this reaction")
                    else:
                        print("Aldehyde is consumed in this reaction")
                    aldehyde_intermediate_detected = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    return aldehyde_intermediate_detected
