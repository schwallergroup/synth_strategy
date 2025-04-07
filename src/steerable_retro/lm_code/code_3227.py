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
    Detects the presence of acetonide-protected diols throughout the synthesis.

    An acetonide protection strategy involves:
    1. The presence of acetonide-protected diols (dioxolane or dioxane rings)
    2. Acetal formation reactions (diol acetalization)
    3. Acetal hydrolysis reactions (to diol)
    """
    acetonide_present = False
    acetonide_formation = False
    acetonide_hydrolysis = False

    def dfs_traverse(node, depth=0):
        nonlocal acetonide_present, acetonide_formation, acetonide_hydrolysis

        # Check molecule nodes for acetonide structures
        if node["type"] == "mol" and "smiles" in node:
            mol_smiles = node["smiles"]

            # Check for dioxolane (5-membered) or dioxane (6-membered) rings
            if checker.check_ring("dioxolane", mol_smiles):
                print(f"Found dioxolane ring (potential acetonide) at depth {depth}: {mol_smiles}")
                acetonide_present = True

            if checker.check_ring("dioxane", mol_smiles):
                print(f"Found dioxane ring (potential acetonide) at depth {depth}: {mol_smiles}")
                acetonide_present = True

        # Check reaction nodes for acetonide formation/hydrolysis
        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rxn_smiles = node["metadata"]["rsmi"]

            # Check for diol acetalization (acetonide formation)
            if checker.check_reaction("Diol acetalization", rxn_smiles):
                print(f"Found diol acetalization reaction at depth {depth}: {rxn_smiles}")
                acetonide_formation = True

            # Check for acetal hydrolysis to diol (acetonide removal)
            if checker.check_reaction("Acetal hydrolysis to diol", rxn_smiles):
                print(f"Found acetal hydrolysis to diol reaction at depth {depth}: {rxn_smiles}")
                acetonide_hydrolysis = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Strategy is present if we have the protected structure or relevant reactions
    strategy_present = acetonide_present or acetonide_formation or acetonide_hydrolysis

    if strategy_present:
        print("Detected acetonide-protected diol strategy:")
        print(f"- Acetonide structures present: {acetonide_present}")
        print(f"- Acetonide formation reactions: {acetonide_formation}")
        print(f"- Acetonide hydrolysis reactions: {acetonide_hydrolysis}")

    return strategy_present
