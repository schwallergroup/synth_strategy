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
    This function detects a synthetic strategy involving C-C bond cleavage,
    specifically methyl group removal.
    """
    has_cc_cleavage = False

    def dfs_traverse(node):
        nonlocal has_cc_cleavage

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1].split(".")

                print(f"Checking reaction: {rsmi}")

                # Check for various C-C bond cleavage reactions
                if checker.check_reaction("Decarboxylation", rsmi):
                    print("Found C-C bond cleavage: Decarboxylation")
                    has_cc_cleavage = True
                elif checker.check_reaction("Oxidation of alkene to carboxylic acid", rsmi):
                    print("Found C-C bond cleavage: Oxidation of alkene to carboxylic acid")
                    has_cc_cleavage = True
                elif checker.check_reaction("Oxidation of ketone to carboxylic acid", rsmi):
                    print("Found C-C bond cleavage: Oxidation of ketone to carboxylic acid")
                    has_cc_cleavage = True
                elif checker.check_reaction(
                    "Ketonization by decarboxylation of carbonic acids", rsmi
                ):
                    print("Found C-C bond cleavage: Ketonization by decarboxylation")
                    has_cc_cleavage = True
                elif checker.check_reaction(
                    "Ketonization by decarboxylation of acid halides", rsmi
                ):
                    print(
                        "Found C-C bond cleavage: Ketonization by decarboxylation of acid halides"
                    )
                    has_cc_cleavage = True
                # Check for oxidative cleavage of alcohols
                elif any(
                    checker.check_fg("Secondary alcohol", r) for r in reactants_smiles
                ) and any(checker.check_fg("Aldehyde", p) for p in product_smiles):
                    print("Found C-C bond cleavage: Secondary alcohol → aldehyde")
                    has_cc_cleavage = True
                elif any(checker.check_fg("Tertiary alcohol", r) for r in reactants_smiles) and any(
                    checker.check_fg("Ketone", p) for p in product_smiles
                ):
                    print("Found C-C bond cleavage: Tertiary alcohol → ketone")
                    has_cc_cleavage = True
                # Check for diol cleavage
                elif (
                    any(checker.check_fg("Primary alcohol", r) for r in reactants_smiles)
                    and any(checker.check_fg("Primary alcohol", r) for r in reactants_smiles)
                    and (
                        any(checker.check_fg("Aldehyde", p) for p in product_smiles)
                        or any(checker.check_fg("Carboxylic acid", p) for p in product_smiles)
                    )
                ):
                    print("Found C-C bond cleavage: Diol cleavage")
                    has_cc_cleavage = True
                # Check for alkene cleavage
                elif any(checker.check_fg("Alkene", r) for r in reactants_smiles) and (
                    any(checker.check_fg("Aldehyde", p) for p in product_smiles)
                    or any(checker.check_fg("Ketone", p) for p in product_smiles)
                    or any(checker.check_fg("Carboxylic acid", p) for p in product_smiles)
                ):
                    print("Found C-C bond cleavage: Alkene cleavage")
                    has_cc_cleavage = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    if has_cc_cleavage:
        print("Detected C-C bond cleavage strategy")
    return has_cc_cleavage
