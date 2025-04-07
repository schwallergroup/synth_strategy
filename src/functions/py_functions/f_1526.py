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
    This function detects if the synthetic route involves a late-stage reduction
    of an aldehyde to an alcohol (typically in the final or penultimate step).

    In retrosynthesis, this means we're looking for an alcohol in the reactants
    that becomes an aldehyde in the products.
    """
    late_stage_reduction_found = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_reduction_found

        # Check reactions at depths 0, 1, or 2 (late stage)
        if node["type"] == "reaction" and depth <= 2:
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]

                # Check if this is an aldehyde reduction reaction
                if checker.check_reaction(
                    "Reduction of aldehydes and ketones to alcohols", rsmi
                ):
                    print(
                        f"Found late-stage aldehyde reduction reaction at depth {depth}"
                    )
                    late_stage_reduction_found = True
                else:
                    # Alternative approach: check for the functional group transformation
                    reactants = rsmi.split(">")[0].split(".")
                    products = rsmi.split(">")[-1].split(".")

                    # In retrosynthesis, we're looking for alcohol in reactants and aldehyde in products
                    alcohol_in_reactants = any(
                        checker.check_fg("Primary alcohol", r)
                        or checker.check_fg("Secondary alcohol", r)
                        or checker.check_fg("Tertiary alcohol", r)
                        for r in reactants
                    )

                    aldehyde_in_products = any(
                        checker.check_fg("Aldehyde", p) for p in products
                    )

                    if alcohol_in_reactants and aldehyde_in_products:
                        print(
                            f"Found late-stage aldehyde reduction (by FG check) at depth {depth}"
                        )
                        late_stage_reduction_found = True

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return late_stage_reduction_found
