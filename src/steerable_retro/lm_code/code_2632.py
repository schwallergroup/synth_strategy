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
    Detects if the synthesis route involves a strategy utilizing multiple electron-withdrawing groups
    (trifluoromethyl and sulfonyl groups) in the final product.
    """
    ewg_strategy_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal ewg_strategy_detected

        if node["type"] == "mol":
            smiles = node["smiles"]

            # Check for trifluoromethyl and sulfonyl groups in the final product (depth 0)
            if depth == 0:
                has_trifluoro = checker.check_fg("Trifluoro group", smiles)
                has_sulfonyl = (
                    checker.check_fg("Sulfonamide", smiles)
                    or checker.check_fg("Sulfone", smiles)
                    or checker.check_fg("Sulfonate", smiles)
                    or checker.check_fg("Sulfonyl halide", smiles)
                )

                if has_trifluoro and has_sulfonyl:
                    print(
                        f"Multiple electron-withdrawing groups strategy detected in final product: {smiles}"
                    )
                    print(f"Contains trifluoromethyl: {has_trifluoro}")
                    print(f"Contains sulfonyl group: {has_sulfonyl}")
                    ewg_strategy_detected = True

        # Traverse children with increased depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root (final product)
    dfs_traverse(route)
    return ewg_strategy_detected
