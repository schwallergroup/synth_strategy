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
    This function detects if the synthetic route involves a late-stage nitro reduction.
    Late stage means it occurs at a low depth in the retrosynthetic tree.

    In retrosynthesis, nitro reduction appears as:
    Product (with amine) -> Reactant (with nitro)
    """
    nitro_reduction_found = False
    nitro_reduction_depth = float("inf")

    def dfs_traverse(node, depth=0):
        nonlocal nitro_reduction_found, nitro_reduction_depth

        if node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if this is a nitro reduction reaction
                if checker.check_reaction("Reduction of nitro groups to amines", rsmi):
                    print(f"Nitro reduction reaction detected at depth {depth}")
                    nitro_reduction_found = True
                    nitro_reduction_depth = min(nitro_reduction_depth, depth)
                else:
                    # Alternative check: look for amine in product and nitro in reactants
                    product_has_amine = (
                        checker.check_fg("Primary amine", product_smiles)
                        or checker.check_fg("Secondary amine", product_smiles)
                        or checker.check_fg("Tertiary amine", product_smiles)
                        or checker.check_fg("Aniline", product_smiles)
                    )

                    if product_has_amine:
                        for reactant_smiles in reactants_smiles:
                            if checker.check_fg("Nitro group", reactant_smiles):
                                # Verify this is likely a nitro reduction by checking atom mapping
                                product_mol = Chem.MolFromSmiles(product_smiles)
                                reactant_mol = Chem.MolFromSmiles(reactant_smiles)

                                if product_mol and reactant_mol:
                                    print(
                                        f"Potential nitro reduction found at depth {depth}"
                                    )
                                    print(f"Product: {product_smiles}")
                                    print(f"Reactant with nitro: {reactant_smiles}")
                                    nitro_reduction_found = True
                                    nitro_reduction_depth = min(
                                        nitro_reduction_depth, depth
                                    )
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Consider it late stage if it's in the first few steps of the synthesis (depth < 3)
    is_late_stage = nitro_reduction_found and nitro_reduction_depth < 3
    print(
        f"Nitro reduction found: {nitro_reduction_found}, at depth: {nitro_reduction_depth}"
    )
    print(f"Is late stage: {is_late_stage}")

    return is_late_stage
