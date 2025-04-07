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
    Checks if the synthesis route contains a sequence of transformations:
    alcohol → halide → azide → amine
    """
    # Track molecules through the sequence
    sequence_molecules = []

    def dfs_traverse(node, depth, path):
        if node["type"] == "mol":
            mol_smiles = node["smiles"]

            # Check if this is a starting material
            if node.get("in_stock", False):
                # Check if starting material has alcohol
                if (
                    checker.check_fg("Primary alcohol", mol_smiles)
                    or checker.check_fg("Secondary alcohol", mol_smiles)
                    or checker.check_fg("Tertiary alcohol", mol_smiles)
                ):
                    path.append(("alcohol", mol_smiles, depth))

            # For intermediate molecules, check for functional groups
            else:
                if (
                    checker.check_fg("Primary halide", mol_smiles)
                    or checker.check_fg("Secondary halide", mol_smiles)
                    or checker.check_fg("Tertiary halide", mol_smiles)
                ):
                    path.append(("halide", mol_smiles, depth))

                if checker.check_fg("Azide", mol_smiles):
                    path.append(("azide", mol_smiles, depth))

                if (
                    checker.check_fg("Primary amine", mol_smiles)
                    or checker.check_fg("Secondary amine", mol_smiles)
                    or checker.check_fg("Tertiary amine", mol_smiles)
                ):
                    path.append(("amine", mol_smiles, depth))

        # Recursively traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1, path)

    # Start traversal
    path = []
    dfs_traverse(route, 0, path)

    # Sort by depth to get chronological order (higher depth = earlier in synthesis)
    path.sort(key=lambda x: x[2], reverse=True)

    # Check if we have the sequence: alcohol → halide → azide → amine
    fg_sequence = [item[0] for item in path]

    # Look for the sequence in order
    try:
        alcohol_idx = fg_sequence.index("alcohol")
        halide_idx = fg_sequence.index("halide")
        azide_idx = fg_sequence.index("azide")
        amine_idx = fg_sequence.index("amine")

        # Check if they appear in the correct order
        if alcohol_idx < halide_idx < azide_idx < amine_idx:
            print(f"Found alcohol→halide→azide→amine sequence: {path}")
            return True
    except ValueError:
        # One of the functional groups not found
        pass

    return False
