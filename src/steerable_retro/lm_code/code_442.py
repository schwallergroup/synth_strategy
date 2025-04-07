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
    Detects if the synthesis route contains a thiazole ring formation reaction.
    """
    has_thiazole_formation = False

    def dfs(node, depth=0):
        nonlocal has_thiazole_formation

        # Check if this is a molecule node with a thiazole ring
        if node["type"] == "mol" and node["smiles"]:
            if checker.check_ring("thiazole", node["smiles"]):
                # Check if this molecule is a product of a reaction
                if not node.get("in_stock", False) and node.get("children", []):
                    for child in node["children"]:
                        if (
                            child["type"] == "reaction"
                            and "metadata" in child
                            and "rsmi" in child["metadata"]
                        ):
                            rxn_smiles = child["metadata"]["rsmi"]
                            product = rxn_smiles.split(">")[-1]
                            reactants = rxn_smiles.split(">")[0].split(".")

                            # Check if thiazole is formed in this reaction (not present in reactants)
                            thiazole_in_reactants = any(
                                checker.check_ring("thiazole", r) for r in reactants
                            )
                            if not thiazole_in_reactants:
                                has_thiazole_formation = True
                                print(f"Found thiazole formation reaction: {rxn_smiles}")

        for child in node.get("children", []):
            dfs(child, depth + 1)

    dfs(route)
    return has_thiazole_formation
