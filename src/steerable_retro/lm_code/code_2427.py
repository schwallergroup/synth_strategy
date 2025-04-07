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
    Detects a synthetic strategy involving early incorporation of a methoxyphenyl group
    that is carried through the synthesis.
    """
    max_depth = 0
    methoxy_depths = []

    def dfs_traverse(node, depth=0):
        nonlocal max_depth

        if depth > max_depth:
            max_depth = depth

        if node["type"] == "mol" and "smiles" in node:
            smiles = node["smiles"]

            # Check for methoxyphenyl pattern using checker
            if checker.check_ring("benzene", smiles) and checker.check_fg("Ether", smiles):
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    # Look for methoxy groups attached to benzene rings
                    pattern = Chem.MolFromSmarts("COc1ccccc1")
                    if mol.HasSubstructMatch(pattern):
                        methoxy_depths.append(depth)
                        print(f"Found methoxyphenyl group at depth {depth} in molecule: {smiles}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    print(f"Max depth: {max_depth}")
    print(f"Methoxyphenyl depths: {methoxy_depths}")

    if not methoxy_depths:
        print("No methoxyphenyl groups found")
        return False

    # Early incorporation means:
    # 1. The group appears in an early stage (high depth value)
    # 2. It persists to the final product (depth 0)
    has_early_methoxy = any(d >= max_depth * 0.6 for d in methoxy_depths)
    has_final_methoxy = 0 in methoxy_depths

    early_incorporation = has_early_methoxy and has_final_methoxy
    print(f"Early incorporation detected: {early_incorporation}")
    return early_incorporation
