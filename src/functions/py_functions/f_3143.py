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
    This function detects if the synthesis route involves esterification of a carboxylic acid.
    """
    esterification_detected = False

    def dfs_traverse(node):
        nonlocal esterification_detected

        if node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if this is an esterification reaction
                if checker.check_reaction("Esterification of Carboxylic Acids", rsmi):
                    print(f"Esterification reaction detected: {rsmi}")
                    esterification_detected = True
                else:
                    # Alternative check: look for carboxylic acid in reactants and ester in product
                    reactants_have_acid = any(
                        checker.check_fg("Carboxylic acid", r) for r in reactants_smiles
                    )
                    product_has_ester = checker.check_fg("Ester", product_smiles)

                    # Check for O-alkylation of carboxylic acids (another esterification method)
                    is_o_alkylation = checker.check_reaction(
                        "O-alkylation of carboxylic acids with diazo compounds", rsmi
                    )

                    if (reactants_have_acid and product_has_ester) or is_o_alkylation:
                        print(f"Carboxylic acid to ester conversion detected: {rsmi}")
                        esterification_detected = True
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return esterification_detected
