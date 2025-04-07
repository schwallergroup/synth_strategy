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


def main(route):
    """
    This function detects a strategy where an aromatic core is maintained throughout
    the synthesis while a side chain is elaborated, eventually leading to cyclization.
    """
    # Track molecules at each step
    molecules_by_depth = {}
    aromatic_core_preserved = False
    side_chain_elaboration = False

    def dfs_traverse(node, depth=0):
        if node["type"] == "mol":
            if "smiles" in node:
                if depth not in molecules_by_depth:
                    molecules_by_depth[depth] = []
                molecules_by_depth[depth].append(node["smiles"])

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Sort depths to analyze molecules in synthesis direction (highest to lowest)
    depths = sorted(molecules_by_depth.keys(), reverse=True)

    if len(depths) < 2:
        return False

    # Check for aromatic core preservation
    try:
        # Get the earliest molecule (highest depth)
        earliest_mol = Chem.MolFromSmiles(molecules_by_depth[depths[0]][0])
        latest_mol = Chem.MolFromSmiles(molecules_by_depth[depths[-1]][0])

        if earliest_mol and latest_mol:
            # Define aromatic pattern
            aromatic_pattern = Chem.MolFromSmarts("c1ccccc1")

            # Check if both earliest and latest molecules have aromatic rings
            if earliest_mol.HasSubstructMatch(
                aromatic_pattern
            ) and latest_mol.HasSubstructMatch(aromatic_pattern):
                aromatic_core_preserved = True
                print("Aromatic core is preserved throughout the synthesis")

            # Check for side chain elaboration
            # Compare the number of non-aromatic atoms
            def count_non_aromatic_atoms(mol):
                count = 0
                for atom in mol.GetAtoms():
                    if not atom.GetIsAromatic():
                        count += 1
                return count

            earliest_non_aromatic = count_non_aromatic_atoms(earliest_mol)
            latest_non_aromatic = count_non_aromatic_atoms(latest_mol)

            if latest_non_aromatic > earliest_non_aromatic:
                side_chain_elaboration = True
                print(
                    f"Side chain elaboration detected: non-aromatic atoms increased from {earliest_non_aromatic} to {latest_non_aromatic}"
                )

    except Exception as e:
        print(f"Error analyzing aromatic core preservation: {e}")

    return aromatic_core_preserved and side_chain_elaboration
