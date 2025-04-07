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
    Detects if a fluorinated aromatic system is maintained throughout the synthesis.
    """
    # Track if we've found a fluorinated aromatic system at each depth
    fluoro_aromatic_at_depth = {}

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            print(f"Processing reaction at depth {depth}")
            # Mark reaction nodes as having the property to maintain continuity
            fluoro_aromatic_at_depth[depth] = True
        elif node["type"] == "mol":
            mol_smiles = node["smiles"]
            print(f"Processing molecule at depth {depth}: {mol_smiles}")

            # Check for fluorinated aromatic
            has_fluoro_aromatic = checker.check_fg("Aromatic halide", mol_smiles)

            # Verify it's specifically a fluorinated aromatic
            if has_fluoro_aromatic:
                mol = Chem.MolFromSmiles(mol_smiles)
                has_fluoro_aromatic = False
                for atom in mol.GetAtoms():
                    if atom.GetSymbol() == "F" and any(
                        neigh.GetIsAromatic() for neigh in atom.GetNeighbors()
                    ):
                        has_fluoro_aromatic = True
                        break

            # Record if this depth has a fluorinated aromatic
            if depth not in fluoro_aromatic_at_depth:
                fluoro_aromatic_at_depth[depth] = has_fluoro_aromatic
            else:
                fluoro_aromatic_at_depth[depth] = (
                    fluoro_aromatic_at_depth[depth] or has_fluoro_aromatic
                )

            print(f"Depth {depth}: {'Has' if has_fluoro_aromatic else 'No'} fluorinated aromatic")

        # Traverse children (reactants)
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if fluorinated aromatic is present at all depths
    if not fluoro_aromatic_at_depth:
        print("No molecules found in route")
        return False

    # Get the maximum depth
    max_depth = max(fluoro_aromatic_at_depth.keys())

    # Check if fluorinated aromatic exists at each depth from 0 to max_depth
    for d in range(max_depth + 1):
        if d not in fluoro_aromatic_at_depth or not fluoro_aromatic_at_depth[d]:
            print(f"No fluorinated aromatic at depth {d}")
            return False

    print(f"Fluorinated aromatic maintained through all {max_depth+1} depths")
    return True
