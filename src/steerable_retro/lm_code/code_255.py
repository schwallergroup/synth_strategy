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
    This function detects if the synthesis maintains a consistent aromatic core
    (benzofuran in this case) throughout the synthesis
    """
    # Track which molecules contain the benzofuran core
    molecules_with_core = 0
    total_molecules = 0

    # Define a benzofuran pattern as backup
    benzofuran_pattern = Chem.MolFromSmarts("c1cc2occc2cc1")

    def dfs_traverse(node):
        nonlocal molecules_with_core, total_molecules

        if node["type"] == "mol" and not node.get("in_stock", False):
            total_molecules += 1
            mol_smiles = node["smiles"]
            print(f"Checking molecule: {mol_smiles}")

            # First try the checker function
            has_benzofuran = checker.check_ring("benzofuran", mol_smiles)

            # If that fails, try RDKit substructure matching
            if not has_benzofuran:
                try:
                    mol = Chem.MolFromSmiles(mol_smiles)
                    if mol and mol.HasSubstructMatch(benzofuran_pattern):
                        has_benzofuran = True
                        print(f"Detected benzofuran via RDKit in: {mol_smiles}")
                except Exception as e:
                    print(f"Error in RDKit processing: {e}")

            # As a last resort, check for common benzofuran SMILES patterns
            if not has_benzofuran:
                if "c1c(oc2c" in mol_smiles or "oc2cc" in mol_smiles:
                    has_benzofuran = True
                    print(f"Detected benzofuran via string matching in: {mol_smiles}")

            if has_benzofuran:
                molecules_with_core += 1
                print(f"Found molecule with benzofuran core: {mol_smiles}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Check if most molecules contain the core
    core_preserved = total_molecules > 0 and molecules_with_core / total_molecules >= 0.8

    print(
        f"Aromatic core preservation detected: {core_preserved} ({molecules_with_core}/{total_molecules} molecules)"
    )
    return core_preserved
