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
    This function detects the presence of a trifluoromethyl group maintained throughout the synthesis.
    """
    # Track CF3 groups at each depth
    cf3_at_depth = {}

    def dfs_traverse(node, current_depth=0):
        if node["type"] == "mol" and "smiles" in node:
            mol_smiles = node["smiles"]

            # Check if molecule contains CF3 group
            if checker.check_fg("Trifluoro group", mol_smiles):
                print(
                    f"Found CF3 group at depth {current_depth} in molecule: {mol_smiles}"
                )
                cf3_at_depth[current_depth] = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, current_depth + 1)

    # Start traversal from root
    dfs_traverse(route)

    # Check if CF3 is present at multiple depths
    if len(cf3_at_depth) >= 2:
        print(
            f"Detected trifluoromethyl group maintained throughout synthesis at depths: {list(cf3_at_depth.keys())}"
        )
        return True

    return False
