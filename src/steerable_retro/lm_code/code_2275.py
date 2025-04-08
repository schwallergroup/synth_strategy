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
    This function detects preservation of a quinolone core structure throughout the synthesis.
    Quinolone core should be present in all intermediate and final products, but not necessarily
    in starting materials.
    """
    all_intermediates_have_quinolone = True

    def dfs_traverse(node, depth=0):
        nonlocal all_intermediates_have_quinolone

        if node["type"] == "mol" and "smiles" in node:
            # Skip starting materials (in_stock)
            if node.get("in_stock", False):
                print(f"Skipping starting material: {node['smiles']}")
                return

            mol_smiles = node["smiles"]

            # Check for quinoline core using the checker function
            has_quinoline = checker.check_ring("quinoline", mol_smiles)
            has_isoquinoline = checker.check_ring("isoquinoline", mol_smiles)

            # Check for carbonyl group attached to the quinoline/isoquinoline
            # (This is what makes it a quinolone)
            mol = Chem.MolFromSmiles(mol_smiles)
            quinolone_pattern = Chem.MolFromSmarts("c1nccc2ccccc12C(=O)")  # Quinoline with carbonyl
            isoquinolone_pattern = Chem.MolFromSmarts(
                "c1ccnc2ccccc12C(=O)"
            )  # Isoquinoline with carbonyl

            has_quinolone = mol.HasSubstructMatch(quinolone_pattern) if mol else False
            has_isoquinolone = mol.HasSubstructMatch(isoquinolone_pattern) if mol else False

            if not (has_quinoline or has_isoquinoline or has_quinolone or has_isoquinolone):
                all_intermediates_have_quinolone = False
                print(f"Found intermediate without quinolone core: {mol_smiles}")
            else:
                print(f"Molecule has quinolone/isoquinolone core: {mol_smiles}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return all_intermediates_have_quinolone
