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
    Detects if the synthesis involves a pyrazolopyrimidine heterocyclic core.
    """
    pyrazolopyrimidine_detected = False

    # SMARTS pattern for pyrazolopyrimidine core
    pyrazolopyrimidine_smarts = "c1nc2c(n1)cnn2"

    def dfs_traverse(node, depth=0):
        nonlocal pyrazolopyrimidine_detected

        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            print(f"Checking molecule at depth {depth}: {mol_smiles}")

            # First check if the molecule contains both pyrazole and pyrimidine rings
            has_pyrazole = checker.check_ring("pyrazole", mol_smiles)
            has_pyrimidine = checker.check_ring("pyrimidine", mol_smiles)

            if has_pyrazole and has_pyrimidine:
                print(
                    f"Found both pyrazole and pyrimidine rings in molecule: {mol_smiles}"
                )

                # Create a molecule object to check for pyrazolopyrimidine core
                mol = Chem.MolFromSmiles(mol_smiles)
                if mol:
                    # Check if the molecule contains a pyrazolopyrimidine core
                    pyrazolopyrimidine_pattern = Chem.MolFromSmarts(
                        pyrazolopyrimidine_smarts
                    )
                    if mol.HasSubstructMatch(pyrazolopyrimidine_pattern):
                        print(f"Pyrazolopyrimidine core detected at depth {depth}")
                        pyrazolopyrimidine_detected = True
                    else:
                        # Alternative approach: check if there's a fused ring system with the right atoms
                        ring_info = mol.GetRingInfo()
                        atom_rings = ring_info.AtomRings()

                        # Look for rings that might be part of a pyrazolopyrimidine
                        for i, ring1 in enumerate(atom_rings):
                            for j, ring2 in enumerate(atom_rings[i + 1 :], i + 1):
                                # Check if rings share atoms (are fused)
                                shared_atoms = set(ring1).intersection(set(ring2))
                                if (
                                    len(shared_atoms) >= 2
                                ):  # Fused rings share at least 2 atoms
                                    # Get atoms in both rings
                                    all_ring_atoms = set(ring1).union(set(ring2))
                                    # Check if this could be a pyrazolopyrimidine
                                    # by counting N atoms (should have at least 3 for pyrazolopyrimidine)
                                    n_count = sum(
                                        1
                                        for atom_idx in all_ring_atoms
                                        if mol.GetAtomWithIdx(atom_idx).GetSymbol()
                                        == "N"
                                    )
                                    if n_count >= 3:
                                        print(
                                            f"Potential pyrazolopyrimidine core detected at depth {depth}"
                                        )
                                        print(
                                            f"Rings share atoms: {shared_atoms}, N count: {n_count}"
                                        )
                                        pyrazolopyrimidine_detected = True

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)
    print(f"Final result: pyrazolopyrimidine_detected = {pyrazolopyrimidine_detected}")
    return pyrazolopyrimidine_detected
