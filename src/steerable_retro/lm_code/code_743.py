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
    This function detects if the synthesis maintains a tetrazole ring throughout.
    """
    # Track tetrazole presence at each depth
    tetrazole_at_depth = {}

    def dfs_traverse(node, depth=0):
        if node["type"] == "mol" and node.get("smiles"):
            mol_smiles = node["smiles"]

            # Check if molecule contains tetrazole using the checker function
            has_tetrazole = checker.check_ring("tetrazole", mol_smiles)

            if has_tetrazole:
                print(f"Tetrazole detected at depth {depth} in molecule: {mol_smiles}")
                tetrazole_at_depth[depth] = mol_smiles

                # Get atom indices of tetrazole rings for more detailed analysis if needed
                tetrazole_indices = checker.get_ring_atom_indices("tetrazole", mol_smiles)
                if tetrazole_indices:
                    print(f"  Tetrazole atom indices: {tetrazole_indices}")

        # Traverse children (reactants in retrosynthetic direction)
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the target molecule
    dfs_traverse(route)

    # Check if tetrazole is present at multiple depths (maintained throughout)
    if len(tetrazole_at_depth) >= 2:
        print(
            f"Tetrazole maintained throughout synthesis at depths: {sorted(tetrazole_at_depth.keys())}"
        )

        # Optional: Print all molecules containing tetrazole
        for depth, smiles in sorted(tetrazole_at_depth.items()):
            print(f"Depth {depth}: {smiles}")

        return True
    else:
        print(
            f"Tetrazole not maintained throughout synthesis. Found at {len(tetrazole_at_depth)} depths."
        )
        return False
