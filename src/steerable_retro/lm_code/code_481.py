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
    This function detects the construction of molecules with 3+ aromatic/heterocyclic rings.
    """
    # Track the maximum number of rings found
    max_ring_count = 0

    def count_rings(mol_smiles):
        if not mol_smiles:
            return 0

        try:
            mol = Chem.MolFromSmiles(mol_smiles)
            if mol is None:
                print(f"Could not parse SMILES: {mol_smiles}")
                return 0

            # Get ring information
            ring_info = mol.GetRingInfo()

            # Count rings - both aromatic and non-aromatic
            total_rings = len(ring_info.AtomRings())

            # Count aromatic rings specifically
            aromatic_rings = 0
            for ring in ring_info.AtomRings():
                if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
                    aromatic_rings += 1

            # Check for heterocyclic rings using the checker function
            heterocyclic_rings = [
                "furan",
                "pyran",
                "pyrrole",
                "pyridine",
                "pyrazole",
                "imidazole",
                "oxazole",
                "thiazole",
                "pyrimidine",
                "pyrazine",
                "pyridazine",
                "triazole",
                "tetrazole",
                "indole",
                "quinoline",
                "isoquinoline",
                "purine",
                "thiophene",
                "benzoxazole",
                "benzothiazole",
                "benzimidazole",
            ]

            hetero_count = 0
            for ring_name in heterocyclic_rings:
                if checker.check_ring(ring_name, mol_smiles):
                    hetero_count += 1
                    print(f"Found {ring_name} ring in molecule")

            # Take the maximum of the different counting methods to avoid undercounting
            # but don't double count by simply adding them
            ring_count = max(total_rings, aromatic_rings, hetero_count)

            print(
                f"Molecule has {ring_count} rings (total: {total_rings}, aromatic: {aromatic_rings}, heterocyclic: {hetero_count})"
            )
            return ring_count

        except Exception as e:
            print(f"Error counting rings: {e}")
            return 0

    def dfs_traverse(node):
        nonlocal max_ring_count

        if node["type"] == "mol" and not node.get("in_stock", False):
            try:
                mol_smiles = node["smiles"]
                ring_count = count_rings(mol_smiles)

                max_ring_count = max(max_ring_count, ring_count)

                if ring_count >= 3:
                    print(f"Found molecule with {ring_count} rings: {mol_smiles}")
            except Exception as e:
                print(f"Error processing molecule node: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    print(f"Maximum ring count found: {max_ring_count}")

    # Return True if we found a molecule with 3+ rings
    return max_ring_count >= 3
