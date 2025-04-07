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
    This function detects chain extension via ethylene oxide addition to form an alcohol.
    """
    ethylene_oxide_used = False

    def dfs_traverse(node):
        nonlocal ethylene_oxide_used

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if any reactant contains oxirane (ethylene oxide ring)
                has_oxirane = any(checker.check_ring("oxirane", r) for r in reactants if r)

                # Check if the product contains a primary alcohol
                has_primary_alcohol = checker.check_fg("Primary alcohol", product)

                # Check if the reaction is a ring opening of epoxide with a nucleophile
                is_epoxide_opening = checker.check_reaction(
                    "Ring opening of epoxide with amine", rsmi
                )

                # If no specific reaction match, check for general epoxide opening pattern
                if not is_epoxide_opening:
                    # Try to create reaction object to check if it's an epoxide opening
                    try:
                        rxn_mol = AllChem.ReactionFromSmarts(rsmi)
                        if (
                            rxn_mol
                            and rxn_mol.GetNumReactantTemplates() > 0
                            and rxn_mol.GetNumProductTemplates() > 0
                        ):
                            # Check if reactants contain oxirane and product doesn't
                            product_has_oxirane = checker.check_ring("oxirane", product)
                            if has_oxirane and not product_has_oxirane and has_primary_alcohol:
                                is_epoxide_opening = True
                    except:
                        pass

                if has_oxirane and has_primary_alcohol and is_epoxide_opening:
                    print(f"Ethylene oxide chain extension detected in reaction: {rsmi}")
                    ethylene_oxide_used = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return ethylene_oxide_used
