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
    Detects if fluorinated aromatic compounds are maintained throughout the synthesis.
    Excludes starting materials from the check.
    """
    all_nodes_have_fluorinated_aromatics = True

    def dfs_traverse(node):
        nonlocal all_nodes_have_fluorinated_aromatics

        if node["type"] == "mol" and node.get("smiles") and not node.get("in_stock", False):
            # Only check non-starting material molecules
            mol_smiles = node["smiles"]

            # Check for fluorinated aromatics using the checker function
            has_aromatic_halide = checker.check_fg("Aromatic halide", mol_smiles)
            has_fluorine = "F" in mol_smiles

            if not (has_aromatic_halide and has_fluorine):
                print(f"Intermediate/product without fluorinated aromatic found: {mol_smiles}")
                all_nodes_have_fluorinated_aromatics = False

        elif node["type"] == "reaction" and node.get("metadata", {}).get("rsmi"):
            # For reaction nodes, check if fluorinated aromatic groups are preserved
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if product has fluorinated aromatics if any reactant has them
            reactants_have_fluorinated = any(
                checker.check_fg("Aromatic halide", r) and "F" in r for r in reactants
            )
            product_has_fluorinated = (
                checker.check_fg("Aromatic halide", product) and "F" in product
            )

            if reactants_have_fluorinated and not product_has_fluorinated:
                print(f"Reaction lost fluorinated aromatic group: {rsmi}")
                all_nodes_have_fluorinated_aromatics = False

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"All nodes have fluorinated aromatics: {all_nodes_have_fluorinated_aromatics}")
    return all_nodes_have_fluorinated_aromatics
