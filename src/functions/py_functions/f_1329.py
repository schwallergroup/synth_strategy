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
    Detects if the synthesis involves multiple chlorinated aromatic rings.
    """
    # Track total unique chlorinated aromatic rings across all molecules
    total_chloro_rings = 0

    def dfs_traverse(node):
        nonlocal total_chloro_rings

        if node["type"] == "mol":
            mol_smiles = node["smiles"]

            # Check if molecule contains aromatic halide functional group
            if checker.check_fg("Aromatic halide", mol_smiles):
                mol = Chem.MolFromSmiles(mol_smiles)
                if mol:
                    # Get ring information
                    ring_info = mol.GetRingInfo()
                    atom_rings = ring_info.AtomRings()

                    # Count chlorinated aromatic rings
                    chloro_ring_count = 0
                    for ring in atom_rings:
                        # Check if ring is aromatic
                        is_aromatic_ring = all(
                            mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring
                        )
                        if not is_aromatic_ring:
                            continue

                        # Check if ring contains a carbon connected to chlorine
                        has_chlorine = False
                        for atom_idx in ring:
                            atom = mol.GetAtomWithIdx(atom_idx)
                            for neighbor in atom.GetNeighbors():
                                if neighbor.GetAtomicNum() == 17:  # Chlorine
                                    has_chlorine = True
                                    break
                            if has_chlorine:
                                break

                        if has_chlorine:
                            chloro_ring_count += 1

                    if chloro_ring_count > 0:
                        print(
                            f"Found molecule with {chloro_ring_count} chlorinated aromatic rings: {mol_smiles}"
                        )
                        total_chloro_rings += chloro_ring_count

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    print(f"Total chlorinated aromatic rings found: {total_chloro_rings}")
    return total_chloro_rings > 1
