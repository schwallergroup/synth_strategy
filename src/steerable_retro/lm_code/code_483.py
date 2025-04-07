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
    Detects a synthesis that maintains an indole core throughout the route.

    This function checks if the indole structure is preserved in the main synthetic pathway,
    excluding starting materials and reagents.
    """
    # Track molecules in the main synthetic pathway
    main_pathway_molecules = []

    def find_main_pathway(node, depth=0):
        """Identify molecules in the main synthetic pathway"""
        if node["type"] == "mol":
            # Skip starting materials (in_stock)
            if not node.get("in_stock", False):
                main_pathway_molecules.append((node["smiles"], depth))

        # Process children
        for child in node.get("children", []):
            find_main_pathway(child, depth + 1)

    # Start traversal to find main pathway
    find_main_pathway(route)

    # Sort by depth (ascending) to get the synthetic pathway from late to early stage
    main_pathway_molecules.sort(key=lambda x: x[1])

    # Check if indole is present in all main pathway molecules
    indole_present_in_all = True

    for mol_smiles, depth in main_pathway_molecules:
        try:
            if not checker.check_ring("indole", mol_smiles):
                indole_present_in_all = False
                print(f"Molecule without indole found in main pathway: {mol_smiles}")
                break
        except Exception as e:
            print(f"Error checking indole in molecule: {e}")
            indole_present_in_all = False
            break

    if indole_present_in_all and main_pathway_molecules:
        print("Indole core is preserved throughout the main synthetic pathway")
    elif not main_pathway_molecules:
        print("No main pathway molecules found to check")
        indole_present_in_all = False

    return indole_present_in_all
