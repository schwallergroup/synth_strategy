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
    Detects a strategy involving the incorporation of a trifluoromethyl-containing fragment
    in a multi-step synthesis with late-stage functionalization.
    """
    # Track if the final product contains a CF3 group
    final_product_has_cf3 = False
    if route["type"] == "mol":
        final_product_has_cf3 = checker.check_fg("Trifluoro group", route["smiles"])
        print(f"Final product has CF3 group: {final_product_has_cf3}")

    # Track if we have a late-stage CF3 incorporation
    has_late_stage_cf3_incorporation = False

    def dfs_traverse(node, depth=0):
        nonlocal has_late_stage_cf3_incorporation

        # Check for CF3 incorporation in reactions
        if (
            node["type"] == "reaction"
            and depth <= 2
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if product has CF3 but not all reactants have it
            product_has_cf3 = checker.check_fg("Trifluoro group", product)
            reactants_with_cf3 = [r for r in reactants if checker.check_fg("Trifluoro group", r)]

            # CF3 incorporation happens when product has CF3 but not all reactants do
            if product_has_cf3 and len(reactants_with_cf3) < len(reactants):
                print(f"Found CF3 incorporation at depth {depth}")

                # Check if this is a late-stage functionalization
                if depth <= 2:
                    print(f"This is a late-stage CF3 incorporation")
                    has_late_stage_cf3_incorporation = True

                    # Check for common CF3 incorporation reactions
                    if checker.check_reaction("Aromatic trifluoromethylation", rsmi):
                        print(f"Detected aromatic trifluoromethylation reaction")
                    elif any("CF3" in r for r in reactants):
                        print(f"Detected CF3-containing fragment incorporation")

        # Recursively process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Strategy is present if final product has CF3 and there was a late-stage CF3 incorporation
    strategy_present = final_product_has_cf3 and has_late_stage_cf3_incorporation
    print(f"Trifluoromethyl-containing fragment strategy detected: {strategy_present}")
    return strategy_present
