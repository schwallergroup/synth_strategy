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
    This function detects if nitrile chemistry is used in the early stages of synthesis (depth 4+).
    """
    early_nitrile_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal early_nitrile_detected

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for nitrile presence in reactants or products
                for compound in reactants:
                    if compound and checker.check_fg("Nitrile", compound):
                        # Check if this is early stage (depth 4+)
                        if depth >= 4:
                            print(
                                f"Detected early-stage nitrile chemistry at depth {depth} in reactant: {compound}"
                            )
                            early_nitrile_detected = True

                # Also check product
                if product and checker.check_fg("Nitrile", product):
                    if depth >= 4:
                        print(
                            f"Detected early-stage nitrile chemistry at depth {depth} in product: {product}"
                        )
                        early_nitrile_detected = True
            except Exception as e:
                print(f"Error processing reaction SMILES at depth {depth}: {e}")

        # Traverse children with incremented depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)
    return early_nitrile_detected
