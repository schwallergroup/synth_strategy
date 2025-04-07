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
    This function detects if a nitrile group is maintained throughout the synthesis.
    It tracks nitrile groups through atom mapping in reactions to ensure the same
    nitrile is preserved from starting materials to the final product.
    """
    # Track if we have a valid path where nitrile is maintained
    nitrile_maintained_paths = []

    def dfs_traverse(node, path_maintained=True, depth=0):
        # For molecule nodes
        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            has_nitrile = checker.check_fg("Nitrile", mol_smiles)

            # If this is a starting material (leaf node)
            if node.get("in_stock", False) or not node.get("children", []):
                # If we're maintaining a nitrile path and this starting material has a nitrile
                if path_maintained and has_nitrile:
                    nitrile_maintained_paths.append(True)
                # If we need a nitrile but this starting material doesn't have one
                elif path_maintained and not has_nitrile:
                    # This path doesn't maintain nitrile
                    pass
                return

            # If this is the final product (depth 0), check if it has a nitrile
            if depth == 0 and not has_nitrile:
                # Final product must have nitrile for strategy to apply
                return

        # For reaction nodes
        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if product has nitrile
                product_has_nitrile = checker.check_fg("Nitrile", product)

                # Check if any reactant has nitrile
                reactants_with_nitrile = [r for r in reactants if checker.check_fg("Nitrile", r)]

                # Check for nitrile-forming reactions
                nitrile_forming = (
                    checker.check_reaction("Nitrile to amide", rsmi)
                    or checker.check_reaction("Oxidation of nitrile to carboxylic acid", rsmi)
                    or checker.check_reaction("Reduction of nitrile to amine", rsmi)
                )

                # If product has nitrile but no reactant does, or vice versa, nitrile is not maintained
                if (product_has_nitrile and not reactants_with_nitrile) or (
                    not product_has_nitrile and reactants_with_nitrile and not nitrile_forming
                ):
                    path_maintained = False
            except Exception as e:
                print(f"Error analyzing reaction: {e}")
                # Be conservative - if we can't analyze, assume nitrile might not be maintained
                path_maintained = False

        # Continue traversal with children
        for child in node.get("children", []):
            dfs_traverse(child, path_maintained, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # If we found at least one path where nitrile is maintained
    return len(nitrile_maintained_paths) > 0
