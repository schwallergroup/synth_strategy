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
    Detects if the synthesis preserves a trifluoromethylpyridine motif throughout.
    """
    molecules_with_motif = 0
    synthesized_molecules = 0

    def dfs_traverse(node):
        nonlocal molecules_with_motif, synthesized_molecules

        if node["type"] == "mol" and node.get("smiles"):
            # Skip in-stock molecules as they are not part of the synthesis
            if node.get("in_stock", False):
                print(f"Skipping in-stock molecule: {node['smiles']}")
                return

            synthesized_molecules += 1
            smiles = node["smiles"]

            try:
                # Check if molecule contains both a pyridine ring and a trifluoro group
                has_pyridine = checker.check_ring("pyridine", smiles)
                has_trifluoro = checker.check_fg("Trifluoro group", smiles)

                print(f"Molecule: {smiles}")
                print(f"Has pyridine: {has_pyridine}, Has trifluoro: {has_trifluoro}")

                # Check if they are connected (trifluoromethylpyridine)
                if has_pyridine and has_trifluoro:
                    mol = Chem.MolFromSmiles(smiles)
                    if not mol:
                        print(f"Could not parse SMILES: {smiles}")
                        return

                    # Get atom indices
                    pyridine_indices = checker.get_ring_atom_indices("pyridine", smiles)
                    trifluoro_indices = checker.get_fg_atom_indices("Trifluoro group", smiles)

                    print(f"Pyridine indices type: {type(pyridine_indices)}")
                    print(f"Pyridine indices: {pyridine_indices}")
                    print(f"Trifluoro indices type: {type(trifluoro_indices)}")
                    print(f"Trifluoro indices: {trifluoro_indices}")

                    # Check if any trifluoro group is directly attached to the pyridine ring
                    is_connected = False

                    # Flatten the indices if they're nested lists/tuples
                    pyridine_atoms = []
                    if pyridine_indices:
                        for match in pyridine_indices:
                            if isinstance(match, tuple) and match:
                                if isinstance(match[0], tuple):
                                    pyridine_atoms.extend(match[0])
                                else:
                                    pyridine_atoms.extend(match)
                            elif isinstance(match, int):
                                pyridine_atoms.append(match)

                    trifluoro_atoms = []
                    if trifluoro_indices:
                        for match in trifluoro_indices:
                            if isinstance(match, tuple) and match:
                                if isinstance(match[0], tuple):
                                    trifluoro_atoms.extend(match[0])
                                else:
                                    trifluoro_atoms.extend(match)
                            elif isinstance(match, int):
                                trifluoro_atoms.append(match)

                    print(f"Flattened pyridine atoms: {pyridine_atoms}")
                    print(f"Flattened trifluoro atoms: {trifluoro_atoms}")

                    # Check for direct connections between any pyridine atom and any trifluoro atom
                    for p_atom in pyridine_atoms:
                        for tf_atom in trifluoro_atoms:
                            try:
                                bond = mol.GetBondBetweenAtoms(p_atom, tf_atom)
                                if bond:
                                    print(
                                        f"Connection found between pyridine atom {p_atom} and trifluoro atom {tf_atom}"
                                    )
                                    is_connected = True
                                    break
                            except Exception as e:
                                print(f"Error checking bond: {e}")
                        if is_connected:
                            break

                    if is_connected:
                        molecules_with_motif += 1
                        print(f"Trifluoromethylpyridine motif detected in: {smiles}")
                    else:
                        # Try an alternative approach - check if any atom in trifluoro is a neighbor of any atom in pyridine
                        for p_atom in pyridine_atoms:
                            p_atom_obj = mol.GetAtomWithIdx(p_atom)
                            for neighbor in p_atom_obj.GetNeighbors():
                                neighbor_idx = neighbor.GetIdx()
                                if neighbor_idx in trifluoro_atoms:
                                    print(
                                        f"Connection found via neighbors: pyridine atom {p_atom} and trifluoro atom {neighbor_idx}"
                                    )
                                    is_connected = True
                                    molecules_with_motif += 1
                                    print(
                                        f"Trifluoromethylpyridine motif detected via neighbors in: {smiles}"
                                    )
                                    break
                            if is_connected:
                                break

                        if not is_connected:
                            print(f"Pyridine and trifluoro present but not connected in: {smiles}")
                else:
                    if not has_pyridine:
                        print(f"No pyridine ring in: {smiles}")
                    if not has_trifluoro:
                        print(f"No trifluoro group in: {smiles}")
            except Exception as e:
                print(f"Error in trifluoromethylpyridine detection: {e}")
                import traceback

                traceback.print_exc()

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    print(f"Synthesized molecules: {synthesized_molecules}")
    print(f"Molecules with trifluoromethylpyridine motif: {molecules_with_motif}")

    # Return True if all synthesized molecules contain the trifluoromethylpyridine motif
    # and we have at least 3 molecules in the route
    return synthesized_molecules >= 3 and molecules_with_motif == synthesized_molecules
