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
    This function detects a strategy involving oxidation of a primary alcohol
    to a carboxylic acid.
    """
    alcohol_to_acid_found = False

    def dfs_traverse(node):
        nonlocal alcohol_to_acid_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this is an oxidation of alcohol to carboxylic acid reaction
            if checker.check_reaction("Oxidation of alcohol to carboxylic acid", rsmi):
                print(f"Detected 'Oxidation of alcohol to carboxylic acid' reaction: {rsmi}")
                alcohol_to_acid_found = True
                return

            # If the specific reaction check fails, check for the functional group transformation
            primary_alcohol_in_reactants = any(
                checker.check_fg("Primary alcohol", r) for r in reactants
            )
            carboxylic_acid_in_product = checker.check_fg("Carboxylic acid", product)
            carboxylic_acid_in_reactants = any(
                checker.check_fg("Carboxylic acid", r) for r in reactants
            )

            # Ensure we have primary alcohol in reactants, carboxylic acid in product,
            # and no carboxylic acid in reactants (to confirm it's a new formation)
            if (
                primary_alcohol_in_reactants
                and carboxylic_acid_in_product
                and not carboxylic_acid_in_reactants
            ):
                print(f"Detected primary alcohol oxidation to carboxylic acid: {rsmi}")
                alcohol_to_acid_found = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return alcohol_to_acid_found
