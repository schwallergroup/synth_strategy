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
    Detects a synthesis strategy where both trifluoromethyl group and
    oxazole heterocycle are preserved throughout the synthesis.
    """
    # Check if target molecule contains the groups
    target_has_cf3 = checker.check_fg("Trifluoro group", route["smiles"])
    target_has_oxazole = checker.check_ring("oxazole", route["smiles"])

    # If target doesn't have both groups, return False
    if not (target_has_cf3 and target_has_oxazole):
        print(
            f"Target molecule doesn't contain both trifluoromethyl and oxazole groups"
        )
        return False

    # Track if each reaction preserves these groups
    all_steps_preserve_cf3 = True
    all_steps_preserve_oxazole = True

    def dfs_traverse(node, depth=0):
        nonlocal all_steps_preserve_cf3, all_steps_preserve_oxazole

        if node["type"] == "mol":
            # Check if molecule contains the groups
            mol_has_cf3 = checker.check_fg("Trifluoro group", node["smiles"])
            mol_has_oxazole = checker.check_ring("oxazole", node["smiles"])

            # Print for debugging
            if mol_has_cf3:
                print(f"Molecule at depth {depth} has trifluoromethyl group")
            if mol_has_oxazole:
                print(f"Molecule at depth {depth} has oxazole ring")

        elif node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check if product contains the groups
            product_has_cf3 = checker.check_fg("Trifluoro group", product_smiles)
            product_has_oxazole = checker.check_ring("oxazole", product_smiles)

            # Check if any reactant contains the groups
            reactants_have_cf3 = any(
                checker.check_fg("Trifluoro group", r) for r in reactants_smiles
            )
            reactants_have_oxazole = any(
                checker.check_ring("oxazole", r) for r in reactants_smiles
            )

            # In retrosynthesis, we check if a group in the product is preserved in at least one reactant
            if product_has_cf3 and not reactants_have_cf3:
                all_steps_preserve_cf3 = False
                print(f"Trifluoromethyl group not preserved at depth {depth}")

            if product_has_oxazole and not reactants_have_oxazole:
                all_steps_preserve_oxazole = False
                print(f"Oxazole heterocycle not preserved at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Return True if both groups are preserved throughout
    return all_steps_preserve_cf3 and all_steps_preserve_oxazole
