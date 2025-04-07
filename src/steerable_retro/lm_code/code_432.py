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
    This function detects a synthetic strategy where a thioether (C-S-C) linkage
    is formed as a key connection between fragments.
    """
    thioether_formed = False

    def dfs_traverse(node):
        nonlocal thioether_formed

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Analyzing reaction: {rsmi}")

                # Check for thiol in reactants
                reactant_has_thiol = any(
                    checker.check_fg("Aliphatic thiol", r) or checker.check_fg("Aromatic thiol", r)
                    for r in reactants
                )
                print(f"Reactant has thiol: {reactant_has_thiol}")

                # Check for alkylating agents in reactants (halides or alcohols)
                reactant_has_alkylating_agent = any(
                    checker.check_fg("Primary halide", r)
                    or checker.check_fg("Secondary halide", r)
                    or checker.check_fg("Tertiary halide", r)
                    or checker.check_fg("Primary alcohol", r)
                    or checker.check_fg("Secondary alcohol", r)
                    or checker.check_fg("Tertiary alcohol", r)
                    for r in reactants
                )
                print(f"Reactant has alkylating agent: {reactant_has_alkylating_agent}")

                # Check for thioether (monosulfide) in product
                product_has_thioether = checker.check_fg("Monosulfide", product)
                print(f"Product has thioether: {product_has_thioether}")

                # Check if reactants don't already have thioether
                reactants_have_thioether = any(
                    checker.check_fg("Monosulfide", r) for r in reactants
                )
                print(f"Reactants have thioether: {reactants_have_thioether}")

                # Check if the reaction is an S-alkylation type or has the characteristics of one
                reaction_is_salkylation = (
                    checker.check_reaction("S-alkylation of thiols", rsmi)
                    or checker.check_reaction("S-alkylation of thiols (ethyl)", rsmi)
                    or checker.check_reaction("S-alkylation of thiols with alcohols", rsmi)
                    or checker.check_reaction("S-alkylation of thiols with alcohols (ethyl)", rsmi)
                    or checker.check_reaction("thioether_nucl_sub", rsmi)
                    or (
                        reactant_has_thiol
                        and reactant_has_alkylating_agent
                        and product_has_thioether
                        and not reactants_have_thioether
                    )
                )
                print(f"Reaction is S-alkylation: {reaction_is_salkylation}")

                # Check if this is a thioether formation (new thioether in product)
                is_thioether_formation = (
                    reactant_has_thiol
                    and reactant_has_alkylating_agent
                    and product_has_thioether
                    and not reactants_have_thioether
                )

                if is_thioether_formation:
                    print(f"Found thioether formation via S-alkylation: {rsmi}")
                    thioether_formed = True

        for child in node.get("children", []):
            dfs_traverse(child)

    # Call dfs_traverse on the root node
    dfs_traverse(route)

    return thioether_formed
