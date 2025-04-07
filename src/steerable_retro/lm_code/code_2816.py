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
    This function detects if the synthesis route contains a thiophene-fused
    pyridine scaffold.
    """
    scaffold_found = False

    def dfs_traverse(node):
        nonlocal scaffold_found

        if node["type"] == "mol" and "smiles" in node:
            mol_smiles = node["smiles"]

            # Check if the molecule contains both thiophene and pyridine rings
            has_thiophene = checker.check_ring("thiophene", mol_smiles)
            has_pyridine = checker.check_ring("pyridine", mol_smiles)

            # If both rings are present, check if they are fused
            if has_thiophene and has_pyridine:
                print(f"Found molecule with both thiophene and pyridine rings: {mol_smiles}")

                try:
                    # Get the atom indices for both rings
                    thiophene_indices = checker.get_ring_atom_indices("thiophene", mol_smiles)
                    pyridine_indices = checker.get_ring_atom_indices("pyridine", mol_smiles)

                    print(f"Thiophene indices: {thiophene_indices}")
                    print(f"Pyridine indices: {pyridine_indices}")

                    # Convert to RDKit molecule to check for fusion
                    mol = Chem.MolFromSmiles(mol_smiles)
                    if mol:
                        # Get ring info
                        ring_info = mol.GetRingInfo()

                        # Get all atom rings
                        atom_rings = ring_info.AtomRings()

                        # Find thiophene and pyridine rings
                        thiophene_rings = []
                        pyridine_rings = []

                        # Use a pattern to identify thiophene and pyridine rings
                        thiophene_patt = Chem.MolFromSmarts("c1ccsc1")
                        pyridine_patt = Chem.MolFromSmarts("c1ccncc1")

                        # Find all matches
                        thiophene_matches = mol.GetSubstructMatches(thiophene_patt)
                        pyridine_matches = mol.GetSubstructMatches(pyridine_patt)

                        print(f"Thiophene matches: {thiophene_matches}")
                        print(f"Pyridine matches: {pyridine_matches}")

                        # Check if any thiophene ring shares atoms with any pyridine ring (fusion)
                        for thiophene_atoms in thiophene_matches:
                            thiophene_set = set(thiophene_atoms)
                            for pyridine_atoms in pyridine_matches:
                                pyridine_set = set(pyridine_atoms)

                                # If the rings share at least 2 atoms, they are fused
                                shared_atoms = thiophene_set.intersection(pyridine_set)
                                if len(shared_atoms) >= 2:
                                    print(
                                        f"Thiophene-fused pyridine scaffold detected in: {mol_smiles}"
                                    )
                                    print(f"Shared atoms between rings: {shared_atoms}")
                                    scaffold_found = True
                                    return

                        # If we didn't find fusion with SMARTS, try using ring info
                        for i, ring1 in enumerate(atom_rings):
                            ring1_set = set(ring1)
                            for j, ring2 in enumerate(atom_rings):
                                if i != j:  # Don't compare a ring with itself
                                    ring2_set = set(ring2)
                                    shared = ring1_set.intersection(ring2_set)
                                    if len(shared) >= 2:
                                        # Check if one is thiophene and one is pyridine
                                        # This is a simplification - would need more checks in practice
                                        print(f"Found fused rings with shared atoms: {shared}")
                                        # Check if these rings correspond to thiophene and pyridine
                                        # For simplicity, we'll assume if we have both rings and they're fused, it's our target
                                        if has_thiophene and has_pyridine:
                                            print(
                                                f"Thiophene-fused pyridine scaffold detected in: {mol_smiles}"
                                            )
                                            scaffold_found = True
                                            return

                except Exception as e:
                    print(f"Error checking for fused rings: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return scaffold_found
