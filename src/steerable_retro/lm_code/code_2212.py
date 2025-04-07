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
    Detects preservation of trifluoromethoxy group throughout the synthesis.
    """
    # Track if the target molecule has a trifluoromethoxy group
    target_has_trifluoromethoxy = False
    # Track if all reactions preserve the trifluoromethoxy group
    all_reactions_preserve_trifluoromethoxy = True

    def dfs_traverse(node, depth=0):
        nonlocal target_has_trifluoromethoxy, all_reactions_preserve_trifluoromethoxy

        if node["type"] == "mol" and "smiles" in node:
            # Check if this molecule has a trifluoromethoxy group
            has_trifluoromethoxy = checker.check_fg("Trifluoro group", node["smiles"])

            # If this is the target molecule (depth 0), record if it has the group
            if depth == 0:
                target_has_trifluoromethoxy = has_trifluoromethoxy
                print(f"Target molecule has trifluoromethoxy: {target_has_trifluoromethoxy}")

        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            # For reaction nodes, check if trifluoromethoxy is preserved
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if product has trifluoromethoxy
                product_has_trifluoromethoxy = checker.check_fg("Trifluoro group", product)

                # Check if any reactant has trifluoromethoxy
                reactants_have_trifluoromethoxy = any(
                    checker.check_fg("Trifluoro group", r) for r in reactants
                )

                # If the product has trifluoromethoxy but no reactant does, it's being introduced
                # If reactants have trifluoromethoxy but product doesn't, it's being removed
                if product_has_trifluoromethoxy and not reactants_have_trifluoromethoxy:
                    print(f"Trifluoromethoxy is being introduced in reaction: {rsmi}")
                elif reactants_have_trifluoromethoxy and not product_has_trifluoromethoxy:
                    print(f"Trifluoromethoxy is being removed in reaction: {rsmi}")
                    all_reactions_preserve_trifluoromethoxy = False
            except Exception as e:
                print(f"Error analyzing reaction: {e}")

        # Traverse children (retrosynthetic direction)
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root (target molecule)
    dfs_traverse(route)

    # The function should return True if:
    # 1. The target molecule has a trifluoromethoxy group
    # 2. All reactions preserve the trifluoromethoxy group (if it exists in the product)
    return target_has_trifluoromethoxy and all_reactions_preserve_trifluoromethoxy
