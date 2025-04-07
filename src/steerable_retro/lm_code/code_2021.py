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
    Checks if the route contains reactions involving morpholine ring manipulation.
    """

    def dfs(node, depth=0):
        if node["type"] == "mol":
            # Check if molecule contains morpholine
            if checker.check_ring("morpholine", node["smiles"]):
                print(f"Found morpholine in molecule: {node['smiles']}")
                return True
        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if reaction involves morpholine
            reactant_has_morpholine = any(checker.check_ring("morpholine", r) for r in reactants)
            product_has_morpholine = checker.check_ring("morpholine", product)

            # If morpholine is created or modified in the reaction
            if (not reactant_has_morpholine and product_has_morpholine) or (
                reactant_has_morpholine and product_has_morpholine
            ):
                print(f"Found morpholine manipulation in reaction: {rsmi}")
                return True

        # Check children
        for child in node.get("children", []):
            if dfs(child, depth + 1):
                return True

        return False

    return dfs(route)
