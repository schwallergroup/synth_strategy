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
    This function detects an azide reduction to amine in the synthesis.
    """
    azide_reduction = False

    def dfs_traverse(node):
        nonlocal azide_reduction

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # First check if this is a known azide reduction reaction
            if checker.check_reaction("Azide to amine reduction (Staudinger)", rsmi):
                print(f"Found azide reduction reaction (Staudinger): {rsmi}")
                azide_reduction = True
            # If not a Staudinger reduction, check for other azide reduction reactions
            elif any(checker.check_fg("Azide", r) for r in reactants) and checker.check_fg(
                "Primary amine", product
            ):
                # Verify that the azide is actually being reduced to an amine
                # by checking that the product doesn't have azide and at least one reactant doesn't have primary amine
                if not checker.check_fg("Azide", product) and any(
                    not checker.check_fg("Primary amine", r) for r in reactants
                ):
                    print(f"Found azide reduction reaction: {rsmi}")
                    azide_reduction = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    return azide_reduction
