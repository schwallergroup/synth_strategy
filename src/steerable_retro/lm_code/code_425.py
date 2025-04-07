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
    This function detects if both an indole heterocycle and a fluorine substituent
    are preserved throughout the synthesis.
    """
    indole_depths = []
    fluorine_depths = []

    def dfs_traverse(node, depth=0):
        nonlocal indole_depths, fluorine_depths

        if node["type"] == "mol" and node.get("smiles"):
            mol_smiles = node["smiles"]
            mol = Chem.MolFromSmiles(mol_smiles)
            if not mol:
                return

            # Check for indole core using the checker function
            if checker.check_ring("indole", mol_smiles):
                indole_depths.append(depth)
                print(f"Indole detected at depth {depth} in molecule: {mol_smiles}")

            # Check specifically for fluorine atoms
            if any(atom.GetSymbol() == "F" for atom in mol.GetAtoms()):
                fluorine_depths.append(depth)
                print(f"Fluorine detected at depth {depth} in molecule: {mol_smiles}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if both motifs are present at multiple depths (preserved throughout)
    indole_preserved = len(indole_depths) >= 2
    fluorine_preserved = len(fluorine_depths) >= 2
    both_preserved = indole_preserved and fluorine_preserved

    print(f"Preserved heterocycle and halogen: {both_preserved}")
    print(f"Indole depths: {indole_depths}, Fluorine depths: {fluorine_depths}")
    return both_preserved
