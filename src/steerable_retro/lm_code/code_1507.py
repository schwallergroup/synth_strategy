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
    This function detects if the synthesis maintains an indole-benzyl scaffold
    throughout the synthesis while modifying other parts of the molecule.
    """
    # Track if indole-benzyl scaffold is present throughout
    scaffold_present_at_all_steps = True
    steps_analyzed = 0

    def dfs_traverse(node):
        nonlocal scaffold_present_at_all_steps, steps_analyzed

        if node["type"] == "mol" and ("in_stock" not in node or not node["in_stock"]):
            mol_smiles = node["smiles"]

            # Check for indole scaffold using checker function
            has_indole = checker.check_ring("indole", mol_smiles)

            # Check for benzyl group using RDKit pattern matching
            # Benzyl group is a phenyl ring connected to a CH2 group
            mol = Chem.MolFromSmiles(mol_smiles)
            benzyl_pattern = Chem.MolFromSmarts("c1ccccc1C")
            has_benzyl = mol.HasSubstructMatch(benzyl_pattern) if mol else False

            print(f"Scaffold check for {mol_smiles}: indole={has_indole}, benzyl={has_benzyl}")

            if not (has_indole and has_benzyl):
                scaffold_present_at_all_steps = False
                print(f"Indole-benzyl scaffold not complete in molecule: {mol_smiles}")
            else:
                print(f"Indole-benzyl scaffold found in molecule: {mol_smiles}")

            steps_analyzed += 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Total steps analyzed: {steps_analyzed}")
    print(f"Scaffold present throughout: {scaffold_present_at_all_steps}")

    # We need to have analyzed at least 3 steps and found the scaffold in all of them
    return scaffold_present_at_all_steps and steps_analyzed >= 3
