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
    This function detects Grignard reactions forming secondary alcohols.
    It looks for reactions where an aldehyde becomes a secondary alcohol.
    """
    grignard_found = False

    def dfs_traverse(node):
        nonlocal grignard_found

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this is a Grignard reaction from aldehyde to alcohol
            if checker.check_reaction("Grignard from aldehyde to alcohol", rsmi):
                print(
                    f"Found Grignard alcohol formation using reaction checker: {rsmi}"
                )
                grignard_found = True
                return

            # Fallback method: check for presence of required functional groups
            has_aldehyde = False
            has_grignard = False
            has_secondary_alcohol = False

            # Check for aldehyde in reactants
            for reactant in reactants:
                if checker.check_fg("Aldehyde", reactant):
                    print(f"Found aldehyde in reactant: {reactant}")
                    has_aldehyde = True

                # Check for Grignard reagent (magnesium halide)
                if checker.check_fg("Magnesium halide", reactant):
                    print(f"Found Grignard reagent in reactant: {reactant}")
                    has_grignard = True

            # Check for secondary alcohol in product
            if checker.check_fg("Secondary alcohol", product):
                print(f"Found secondary alcohol in product: {product}")
                has_secondary_alcohol = True

            # A Grignard reaction requires both an aldehyde and a Grignard reagent
            if has_aldehyde and has_grignard and has_secondary_alcohol:
                print(
                    f"Found Grignard alcohol formation using functional group checks: {rsmi}"
                )
                grignard_found = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(f"Grignard alcohol formation found: {grignard_found}")
    return grignard_found
