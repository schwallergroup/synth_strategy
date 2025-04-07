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
    This function detects if the synthesis involves an indole core structure.
    """
    indole_found = False

    def dfs_traverse(node):
        nonlocal indole_found

        if indole_found:
            return  # Early return if indole already found

        if node["type"] == "mol" and node.get("in_stock", False) == False:
            smiles = node["smiles"]

            # Check for indole core using the checker function
            try:
                if checker.check_ring("indole", smiles):
                    print(f"Found molecule with indole core: {smiles}")
                    indole_found = True
            except Exception as e:
                print(f"Error checking for indole in molecule {smiles}: {e}")

        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            # Check reaction components for indole
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check product for indole
                if checker.check_ring("indole", product):
                    print(f"Found product with indole core in reaction: {product}")
                    indole_found = True

                # Check reactants for indole
                for reactant in reactants:
                    if checker.check_ring("indole", reactant):
                        print(f"Found reactant with indole core in reaction: {reactant}")
                        indole_found = True
            except Exception as e:
                print(f"Error checking reaction components for indole: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return indole_found
