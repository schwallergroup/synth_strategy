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
    Detects preservation of bis-pyrazole scaffold throughout the synthesis
    """
    # Track if pyrazole rings are present at each step
    all_steps_have_pyrazoles = True
    step_count = 0

    def dfs_traverse(node):
        nonlocal all_steps_have_pyrazoles, step_count

        if node["type"] == "mol" and "smiles" in node:
            mol_smiles = node["smiles"]

            # Skip checking reagents/catalysts that are typically not part of the main scaffold
            if mol_smiles in ["O=S(Cl)Cl", "N", "O=S(=O)(Cl)Cl"]:
                print(f"Skipping reagent/catalyst: {mol_smiles}")
                return

            # Check for pyrazole rings using the checker function
            if checker.check_ring("pyrazole", mol_smiles):
                # Count pyrazole rings
                pyrazole_indices = checker.get_ring_atom_indices("pyrazole", mol_smiles)
                print(f"Found {len(pyrazole_indices)} pyrazole rings in molecule: {mol_smiles}")

                if len(pyrazole_indices) < 2:  # We expect at least 2 pyrazole rings
                    all_steps_have_pyrazoles = False
                    print(
                        f"Found only {len(pyrazole_indices)} pyrazole rings in molecule: {mol_smiles}"
                    )
            else:
                # Try to check if the molecule contains connected pyrazoles that might not be detected individually
                mol = Chem.MolFromSmiles(mol_smiles)
                if "cnn1-c1" in mol_smiles or "nn(C)c" in mol_smiles:
                    print(f"Detected connected pyrazole structure in: {mol_smiles}")
                    # Count how many pyrazole-like structures we have
                    pyrazole_count = mol_smiles.count("cnn") + mol_smiles.count("nn(")
                    if pyrazole_count >= 2:
                        print(f"Found approximately {pyrazole_count} pyrazole-like structures")
                    else:
                        all_steps_have_pyrazoles = False
                        print(f"Insufficient pyrazole-like structures: {pyrazole_count}")
                else:
                    all_steps_have_pyrazoles = False
                    print(f"No pyrazole rings found in molecule: {mol_smiles}")

            step_count += 1

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Only return True if we've checked at least one molecule and all have pyrazoles
    return all_steps_have_pyrazoles and step_count > 0
