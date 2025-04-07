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
    Detects if the synthesis route involves a dibrominated aromatic intermediate.
    """
    dibrominated_found = False

    def dfs_traverse(node):
        nonlocal dibrominated_found

        if node["type"] == "mol":
            smiles = node["smiles"]
            # First check if the molecule contains aromatic halides
            if checker.check_fg("Aromatic halide", smiles):
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    # Count bromine atoms attached to aromatic carbons
                    bromine_count = 0
                    for atom in mol.GetAtoms():
                        # Check if atom is bromine
                        if atom.GetSymbol() == "Br":
                            # Get the atom it's connected to
                            neighbors = atom.GetNeighbors()
                            if neighbors and len(neighbors) > 0:
                                # Check if the neighbor is aromatic
                                if neighbors[0].GetIsAromatic():
                                    bromine_count += 1

                    # If there are at least 2 bromine atoms on aromatic rings, it's dibrominated
                    if bromine_count >= 2:
                        print(f"Dibrominated aromatic intermediate found: {smiles}")
                        print(f"Number of bromine atoms on aromatic rings: {bromine_count}")
                        dibrominated_found = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    print(f"Dibrominated aromatic intermediate strategy: {dibrominated_found}")
    return dibrominated_found
