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
    This function detects the activation of an alcohol by conversion to a mesylate
    as a leaving group preparation strategy.
    """

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if any reactant contains an alcohol group
            reactant_has_alcohol = any(
                checker.check_fg("Primary alcohol", r)
                or checker.check_fg("Secondary alcohol", r)
                or checker.check_fg("Tertiary alcohol", r)
                or checker.check_fg("Aromatic alcohol", r)
                for r in reactants
                if r
            )

            # Check if the product contains a mesylate group
            product_has_mesylate = checker.check_fg("Mesylate", product) if product else False

            # Check for methanesulfonyl chloride or similar reagent
            has_mesylating_agent = any("S(=O)(=O)Cl" in r for r in reactants if r)

            # Check if this is a sulfonic ester formation reaction
            is_sulfonic_ester_formation = checker.check_reaction(
                "Formation of Sulfonic Esters", rsmi
            ) or checker.check_reaction(
                "Formation of Sulfonic Esters on TMS protected alcohol", rsmi
            )

            print(
                f"Depth: {depth}, Alcohol: {reactant_has_alcohol}, Mesylate: {product_has_mesylate}, "
                f"Mesylating agent: {has_mesylating_agent}, Sulfonic ester formation: {is_sulfonic_ester_formation}"
            )

            if reactant_has_alcohol and product_has_mesylate and is_sulfonic_ester_formation:
                print(f"Found alcohol activation to mesylate at depth {depth}")
                return True

            # Alternative check with less strict requirements
            if reactant_has_alcohol and product_has_mesylate:
                print(
                    f"Found potential alcohol activation to mesylate at depth {depth} (reaction type not confirmed)"
                )
                return True

        for child in node.get("children", []):
            if dfs_traverse(child, depth + 1):
                return True

        return False

    result = dfs_traverse(route)
    return result
