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
    Detects phthalimide protection strategy in the synthetic route.
    """
    phthalimide_found = False

    def dfs(node, depth=0):
        nonlocal phthalimide_found

        if phthalimide_found:
            return

        if node["type"] == "mol":
            # Check if molecule contains phthalimide group
            if checker.check_fg("Unsubstituted dicarboximide", node["smiles"]) or checker.check_fg(
                "Substituted dicarboximide", node["smiles"]
            ):
                print(f"Found phthalimide group in molecule: {node['smiles']}")
                phthalimide_found = True

        elif node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for phthalimide protection reaction
                if checker.check_reaction("Phthalic anhydride to phthalimide", rsmi):
                    print(f"Found phthalimide protection reaction: {rsmi}")
                    phthalimide_found = True

                # Also check for phthalimide deprotection
                if checker.check_reaction("Phthalimide deprotection", rsmi):
                    print(f"Found phthalimide deprotection reaction: {rsmi}")
                    phthalimide_found = True
            except Exception as e:
                print(f"Error in phthalimide check: {e}")

        # Continue DFS traversal
        for child in node.get("children", []):
            dfs(child, depth + 1)

    dfs(route)
    print(f"Phthalimide protection strategy detected: {phthalimide_found}")
    return phthalimide_found
