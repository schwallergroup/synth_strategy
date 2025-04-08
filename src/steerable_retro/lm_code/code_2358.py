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
    Detects thioether fragment coupling strategy in the synthesis route.
    """
    thioether_formation_found = False

    def dfs(node, depth=0):
        nonlocal thioether_formation_found

        if node["type"] == "reaction":
            try:
                rxn_smiles = node.get("metadata", {}).get("rsmi", "")
                if not rxn_smiles:
                    return

                # Check for thioether formation reactions
                if (
                    checker.check_reaction("S-alkylation of thiols", rxn_smiles)
                    or checker.check_reaction("S-alkylation of thiols (ethyl)", rxn_smiles)
                    or checker.check_reaction("S-alkylation of thiols with alcohols", rxn_smiles)
                    or checker.check_reaction(
                        "S-alkylation of thiols with alcohols (ethyl)", rxn_smiles
                    )
                    or checker.check_reaction("thioether_nucl_sub", rxn_smiles)
                ):

                    # Verify that a thioether is actually formed
                    reactants = rxn_smiles.split(">")[0].split(".")
                    product = rxn_smiles.split(">")[-1]

                    # Check if any reactant has a thiol group
                    has_thiol = any(
                        checker.check_fg("Aliphatic thiol", r)
                        or checker.check_fg("Aromatic thiol", r)
                        for r in reactants
                    )

                    # Check if product has a monosulfide (thioether) group
                    has_thioether = checker.check_fg("Monosulfide", product)

                    if has_thiol and has_thioether:
                        thioether_formation_found = True
                        print(f"Found thioether formation reaction: {rxn_smiles}")

                # Additional check for thioether formation through other reactions
                if not thioether_formation_found:
                    reactants = rxn_smiles.split(">")[0].split(".")
                    product = rxn_smiles.split(">")[-1]

                    # Check if any reactant has a thiol group and product has a thioether
                    has_thiol = any(
                        checker.check_fg("Aliphatic thiol", r)
                        or checker.check_fg("Aromatic thiol", r)
                        for r in reactants
                    )
                    has_thioether = checker.check_fg("Monosulfide", product)

                    if (
                        has_thiol
                        and has_thioether
                        and not any(checker.check_fg("Monosulfide", r) for r in reactants)
                    ):
                        thioether_formation_found = True
                        print(
                            f"Found thioether formation through alternative reaction: {rxn_smiles}"
                        )

            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Recursively check children
        for child in node.get("children", []):
            dfs(child, depth + 1)

    dfs(route)
    print(f"Thioether formation found: {thioether_formation_found}")
    return thioether_formation_found
