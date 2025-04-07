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
    Detects if fluorophenyl groups are maintained throughout the synthesis.
    """
    depths_with_fluorophenyl = set()
    mol_depths = set()
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal depths_with_fluorophenyl, mol_depths, max_depth

        if node["type"] == "mol":
            mol_depths.add(depth)
            mol_smiles = node["smiles"]

            # Check if the molecule contains a fluorophenyl group
            # First check for aromatic halide, then confirm it's specifically fluorine
            if checker.check_fg("Aromatic halide", mol_smiles) and "F" in mol_smiles:
                # Verify it's actually a fluorophenyl by checking the molecule
                mol = Chem.MolFromSmiles(mol_smiles)
                if mol:
                    # Look for aromatic carbon connected to fluorine
                    fluorophenyl_pattern = Chem.MolFromSmarts("c-F")
                    if mol.HasSubstructMatch(fluorophenyl_pattern):
                        depths_with_fluorophenyl.add(depth)
                        print(
                            f"Fluorophenyl group detected at depth {depth} in molecule: {mol_smiles}"
                        )

        max_depth = max(max_depth, depth)

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)

    # Check if fluorophenyl groups are present at all molecule depths
    if not mol_depths:
        print("No molecule nodes found in the route")
        return False

    all_mol_depths_with_fluorophenyl = all(
        depth in depths_with_fluorophenyl for depth in mol_depths
    )

    if all_mol_depths_with_fluorophenyl:
        print("Fluorophenyl groups maintained throughout the synthesis")
    else:
        missing_depths = [d for d in mol_depths if d not in depths_with_fluorophenyl]
        print(f"Fluorophenyl groups missing at depths: {missing_depths}")

    return all_mol_depths_with_fluorophenyl
