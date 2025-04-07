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
    This function detects if the synthetic route produces a compound containing
    a diarylamine motif.
    """
    contains_diarylamine = False

    def is_diarylamine(smiles):
        """Check if a molecule contains a diarylamine (secondary amine connecting two aromatic rings)"""
        # First check if it contains a secondary amine
        if not checker.check_fg("Secondary amine", smiles):
            return False

        # Then verify that the secondary amine connects two aromatic rings
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return False

        # Find all secondary amine nitrogens
        sec_amine_pattern = Chem.MolFromSmarts("[NX3;H1;!$(NC=O)]")
        matches = mol.GetSubstructMatches(sec_amine_pattern)

        for match in matches:
            n_idx = match[0]
            n_atom = mol.GetAtomWithIdx(n_idx)

            # Get the atoms connected to the nitrogen
            neighbors = [a for a in n_atom.GetNeighbors()]
            if len(neighbors) != 2:
                continue

            # Check if both neighbors are aromatic carbons
            if (
                neighbors[0].GetIsAromatic()
                and neighbors[0].GetSymbol() == "C"
                and neighbors[1].GetIsAromatic()
                and neighbors[1].GetSymbol() == "C"
            ):
                print(f"Found diarylamine in {smiles}")
                return True

        return False

    def dfs_traverse(node, depth=0):
        nonlocal contains_diarylamine

        if node["type"] == "mol":
            # Check if this molecule contains a diarylamine
            smiles = node.get("smiles", "")
            print(f"Checking molecule at depth {depth}: {smiles}")

            if is_diarylamine(smiles):
                contains_diarylamine = True
                print(f"Diarylamine motif detected in molecule: {smiles}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    print("Starting traversal of synthetic route")
    dfs_traverse(route)
    print(f"Diarylamine detection result: {contains_diarylamine}")

    return contains_diarylamine
