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
    Detects if the synthesis route contains a nitro reduction to amine reaction.
    """
    has_nitro_reduction = False

    def dfs(node, depth=0):
        nonlocal has_nitro_reduction

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rxn_smiles = node["metadata"]["rsmi"]

            # Check if this is a nitro reduction reaction
            if checker.check_reaction("Reduction of nitro groups to amines", rxn_smiles):
                has_nitro_reduction = True
                print(f"Found nitro reduction reaction: {rxn_smiles}")
            else:
                # Alternative check: look for nitro group in reactants and amine in products
                reactants = rxn_smiles.split(">")[0].split(".")
                product = rxn_smiles.split(">")[-1]

                nitro_in_reactants = any(checker.check_fg("Nitro group", r) for r in reactants)
                amine_in_product = (
                    checker.check_fg("Primary amine", product)
                    or checker.check_fg("Secondary amine", product)
                    or checker.check_fg("Tertiary amine", product)
                )

                if (
                    nitro_in_reactants
                    and amine_in_product
                    and not checker.check_fg("Nitro group", product)
                ):
                    has_nitro_reduction = True
                    print(f"Found nitro reduction reaction (by FG analysis): {rxn_smiles}")

        for child in node.get("children", []):
            dfs(child, depth + 1)

    dfs(route)
    return has_nitro_reduction
