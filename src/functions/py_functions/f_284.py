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
    This function detects nitro reduction as a key functional group interconversion.
    """
    found_nitro_reduction = False

    def dfs_traverse(node, depth=0):
        nonlocal found_nitro_reduction

        if node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if this is a nitro reduction reaction
                if checker.check_reaction("Reduction of nitro groups to amines", rsmi):
                    found_nitro_reduction = True
                    print(
                        f"Nitro reduction detected at depth {depth} with reaction SMILES: {rsmi}"
                    )
                    return

                # Alternative approach: check for nitro group in reactants and amine in product
                for reactant in reactants_smiles:
                    if checker.check_fg("Nitro group", reactant):
                        # If reactant has nitro group, check if product has amine
                        if checker.check_fg("Primary amine", product_smiles):
                            # Verify it's a reduction reaction (could be more specific)
                            found_nitro_reduction = True
                            print(
                                f"Nitro reduction detected at depth {depth} - Nitro in reactant and amine in product"
                            )
                            return
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return found_nitro_reduction
