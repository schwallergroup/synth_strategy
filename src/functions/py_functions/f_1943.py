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
    This function detects if a 1,3-dichloroaromatic motif is preserved throughout the synthesis.
    """
    all_intermediates_have_dichloro = True

    def dfs_traverse(node, depth=0):
        nonlocal all_intermediates_have_dichloro

        if node["type"] == "mol" and node.get("smiles"):
            # Skip checking starting materials (in_stock)
            if not node.get("in_stock", False):
                # Check if the molecule has a 1,3-dichloroaromatic pattern
                mol_smiles = node["smiles"]
                # Use a more accurate pattern for 1,3-dichloroaromatic motif
                has_dichloro = False

                # Check for 1,3-dichlorobenzene pattern
                if checker.check_fg("Aromatic halide", mol_smiles):
                    mol = Chem.MolFromSmiles(mol_smiles)
                    if mol:
                        # Check specifically for 1,3-dichloro pattern on aromatic ring
                        dichloro_pattern = Chem.MolFromSmarts("c1(Cl)cc(Cl)ccc1")
                        if mol.HasSubstructMatch(dichloro_pattern):
                            has_dichloro = True

                if not has_dichloro:
                    all_intermediates_have_dichloro = False
                    print(
                        f"Intermediate/product without 1,3-dichloro pattern found: {mol_smiles}"
                    )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return all_intermediates_have_dichloro
