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
    This function detects a functional group interconversion from primary amine to nitrile.
    In retrosynthesis, this means detecting nitrile to primary amine transformation.
    """
    # Track if we've seen the transformation
    amine_to_nitrile_found = False

    def dfs_traverse(node, depth=0):
        nonlocal amine_to_nitrile_found

        if node["type"] == "reaction":
            try:
                # Extract reactants and product (retrosynthetically)
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # In retrosynthesis, we're moving from target to starting materials
                # So we check for nitrile in reactants and primary amine in product
                reactants_have_nitrile = any(
                    checker.check_fg("Nitrile", r) for r in reactants_smiles
                )

                product_has_amine = checker.check_fg("Primary amine", product_smiles)

                # If reactants have nitrile and product has primary amine, transformation occurred
                if reactants_have_nitrile and product_has_amine:
                    print(
                        f"Potential nitrile to amine transformation detected at depth {depth}"
                    )

                    # Check for known reactions that convert nitriles to amines
                    if checker.check_reaction("Reduction of nitrile to amine", rsmi):
                        print(f"Nitrile to amine reduction detected at depth {depth}")
                        amine_to_nitrile_found = True
                    else:
                        # Even if we don't recognize the specific reaction, if we have confirmed
                        # the functional group change, we'll count it
                        print(
                            f"Unclassified nitrile to amine transformation detected at depth {depth}"
                        )
                        amine_to_nitrile_found = True
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return amine_to_nitrile_found
