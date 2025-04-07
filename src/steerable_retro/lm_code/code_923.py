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
    This function detects if the synthesis involves a C-S bond formation
    via coupling of an aryl halide with a thiol.
    """
    found_coupling = False

    def dfs_traverse(node, depth=0):
        nonlocal found_coupling

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0]
            product_smiles = rsmi.split(">")[-1]

            try:
                # First check if the reaction is directly an Ullmann-Goldberg thiol substitution
                if checker.check_reaction("Ullmann-Goldberg Substitution thiol", rsmi):
                    found_coupling = True
                    print(f"Found Ullmann-Goldberg thiol substitution at depth {depth}")
                    return

                # If not, check for the components manually
                reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles.split(".")]
                product = Chem.MolFromSmiles(product_smiles)

                # Check for aryl halide in reactants
                has_aryl_halide = any(
                    r and checker.check_fg("Aromatic halide", Chem.MolToSmiles(r))
                    for r in reactants
                    if r
                )

                # Check for thiol in reactants (both aromatic and aliphatic)
                has_thiol = any(
                    r
                    and (
                        checker.check_fg("Aromatic thiol", Chem.MolToSmiles(r))
                        or checker.check_fg("Aliphatic thiol", Chem.MolToSmiles(r))
                    )
                    for r in reactants
                    if r
                )

                # Check for thioether (monosulfide) in product
                has_thioether_product = product and checker.check_fg(
                    "Monosulfide", Chem.MolToSmiles(product)
                )

                # Check if all conditions are met
                if has_aryl_halide and has_thiol and has_thioether_product:
                    # Additional check to ensure C-S bond formation
                    found_coupling = True
                    print(f"Found aryl halide-thiol coupling at depth {depth}")
            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Aryl halide-thiol coupling detection result: {found_coupling}")
    return found_coupling
