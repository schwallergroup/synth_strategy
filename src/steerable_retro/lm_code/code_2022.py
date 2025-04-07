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


def main(route, min_count=2):
    """
    Checks if the route contains multiple SNAr (nucleophilic aromatic substitution) reactions.
    """
    snar_count = 0

    def dfs(node, depth=0):
        nonlocal snar_count

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]

            # Check for SNAr reactions
            # Look for nucleophilic aromatic substitution patterns
            if (
                checker.check_reaction("Ullmann-Goldberg Substitution amine", rsmi)
                or checker.check_reaction("Ullmann-Goldberg Substitution thiol", rsmi)
                or checker.check_reaction("Ullmann-Goldberg Substitution aryl alcohol", rsmi)
                or checker.check_reaction("heteroaromatic_nuc_sub", rsmi)
                or checker.check_reaction("nucl_sub_aromatic_ortho_nitro", rsmi)
                or checker.check_reaction("nucl_sub_aromatic_para_nitro", rsmi)
            ):
                snar_count += 1
                print(f"Found SNAr reaction ({snar_count}): {rsmi}")
                if snar_count >= min_count:
                    return True

        # Check children
        for child in node.get("children", []):
            if dfs(child, depth + 1):
                return True

        return False

    dfs(route)
    return snar_count >= min_count
