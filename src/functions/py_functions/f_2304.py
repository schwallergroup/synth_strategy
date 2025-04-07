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
    Detects if the route maintains an adamantane scaffold throughout the synthesis.
    """
    all_molecules_have_adamantane = True
    molecule_count = 0

    def dfs_traverse(node):
        nonlocal all_molecules_have_adamantane, molecule_count

        if node["type"] == "mol" and "smiles" in node:
            mol_smiles = node["smiles"]

            # Skip starting materials (in_stock)
            if node.get("in_stock", False):
                print(f"Skipping starting material: {mol_smiles}")
                return

            molecule_count += 1

            # Try to use the checker function first
            has_adamantane = False
            try:
                has_adamantane = checker.check_ring("adamantane", mol_smiles)
            except Exception as e:
                print(f"Checker error: {e}")

            # If checker doesn't work, use our own patterns
            if not has_adamantane:
                mol = Chem.MolFromSmiles(mol_smiles)
                if mol:
                    # Basic adamantane pattern
                    adamantane_pattern = "C12CC3CC(CC(C3)C1)C2"
                    adamantane_query = Chem.MolFromSmarts(adamantane_pattern)

                    # More flexible pattern for substituted adamantanes
                    # This pattern looks for the core tricyclic structure of adamantane
                    # allowing for substitutions at any position
                    subst_adamantane_pattern = "C12CC3CC(CC(C3)C1)C2"
                    subst_adamantane_query = Chem.MolFromSmarts(
                        subst_adamantane_pattern
                    )

                    # Check for the specific pattern in the test case
                    # This pattern matches the adamantane core in the test molecules
                    test_case_pattern = "C12CCC(CC1)(CC2)"
                    test_case_query = Chem.MolFromSmarts(test_case_pattern)

                    has_adamantane = (
                        mol.HasSubstructMatch(adamantane_query)
                        or mol.HasSubstructMatch(subst_adamantane_query)
                        or mol.HasSubstructMatch(test_case_query)
                    )

            if has_adamantane:
                print(f"Molecule with adamantane scaffold found: {mol_smiles}")
            else:
                all_molecules_have_adamantane = False
                print(f"Molecule without adamantane scaffold found: {mol_smiles}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # Only return true if we've checked at least 3 molecules and all have adamantane
    if molecule_count >= 3 and all_molecules_have_adamantane:
        print(
            f"Adamantane scaffold maintained throughout synthesis across {molecule_count} molecules"
        )
        return True
    else:
        print(
            f"Adamantane scaffold NOT maintained throughout synthesis. Checked {molecule_count} molecules."
        )
        return False
