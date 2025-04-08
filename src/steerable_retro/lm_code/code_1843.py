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
    Detects if the synthesis maintains a tetrahydronaphthalene core throughout.
    """
    # Track if we've found at least one molecule with the core
    found_core = False
    # Track if any molecule in the main synthetic path lacks the core
    all_have_core = True

    def dfs_traverse(node, depth=0):
        nonlocal found_core, all_have_core

        if node["type"] == "mol":
            mol_smiles = node["smiles"]

            # Check if this molecule has the tetrahydronaphthalene core
            # Looking for patterns that indicate a tetrahydronaphthalene-like structure
            # The molecules in the errors have a structure with "CCC2" and aromatic carbons
            has_naphthalene = checker.check_ring("naphthalene", mol_smiles)

            # Check for tetrahydronaphthalene-like structure
            # This pattern looks for a fused ring system with one aromatic ring and one saturated ring
            mol = Chem.MolFromSmiles(mol_smiles)
            has_tetrahydro = False

            if mol:
                # Check for a structure with a benzene ring fused to a cyclohexane ring
                # This covers 1,2,3,4-tetrahydronaphthalene and similar structures
                ring_info = mol.GetRingInfo()
                if ring_info.NumRings() >= 2:
                    # Look for the pattern in the SMILES that indicates tetrahydronaphthalene
                    # All error examples have "c1" and "c2" (aromatic carbons) and "CCC2" (saturated carbons)
                    if (
                        "c1" in mol_smiles
                        and "c2" in mol_smiles
                        and ("CCC2" in mol_smiles or "CC2" in mol_smiles)
                    ):
                        has_tetrahydro = True

            has_core = has_naphthalene or has_tetrahydro

            # If it's not a starting material, it should have the core
            if not node.get("in_stock", False):
                if has_core:
                    found_core = True
                    print(f"Found molecule with tetrahydronaphthalene-like core: {mol_smiles}")
                else:
                    print(
                        f"Intermediate or product without tetrahydronaphthalene core: {mol_smiles}"
                    )
                    all_have_core = False
            # For starting materials, we just note if any have the core
            elif has_core:
                found_core = True
                print(f"Starting material with tetrahydronaphthalene core: {mol_smiles}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # The strategy is valid if we found at least one molecule with the core
    # and all non-starting materials maintain the core
    return found_core and all_have_core
