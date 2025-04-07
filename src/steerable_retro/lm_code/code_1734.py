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
    Detects if the synthesis route preserves specific functional groups (nitrile and trifluoromethyl)
    throughout the synthesis.
    """
    # Track if we've found both functional groups in the target molecule
    target_has_nitrile = False
    target_has_trifluoromethyl = False

    # Check if target molecule has both functional groups
    if route["type"] == "mol":
        target_mol_smiles = route["smiles"]
        target_has_nitrile = checker.check_fg("Nitrile", target_mol_smiles)
        target_has_trifluoromethyl = checker.check_fg("Trifluoro group", target_mol_smiles)

    # If target doesn't have both functional groups, strategy doesn't apply
    if not (target_has_nitrile and target_has_trifluoromethyl):
        print("Target molecule doesn't have both nitrile and trifluoromethyl groups")
        return False

    # Track preservation through reactions
    preservation = True

    def check_reaction_preservation(node):
        nonlocal preservation

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check if product has the functional groups
            product_has_nitrile = checker.check_fg("Nitrile", product_smiles)
            product_has_trifluoromethyl = checker.check_fg("Trifluoro group", product_smiles)

            # If product has the functional groups, at least one reactant should have them too
            if product_has_nitrile:
                reactants_with_nitrile = any(
                    checker.check_fg("Nitrile", r) for r in reactants_smiles
                )
                if not reactants_with_nitrile:
                    print(f"Nitrile group created in reaction: {rsmi}")
                    preservation = False

            if product_has_trifluoromethyl:
                reactants_with_trifluoromethyl = any(
                    checker.check_fg("Trifluoro group", r) for r in reactants_smiles
                )
                if not reactants_with_trifluoromethyl:
                    print(f"Trifluoromethyl group created in reaction: {rsmi}")
                    preservation = False

        # Traverse children
        for child in node.get("children", []):
            check_reaction_preservation(child)

    # Start traversal from the root
    check_reaction_preservation(route)

    if preservation:
        print("Both nitrile and trifluoromethyl groups are preserved throughout the synthesis")

    return preservation
