#!/bin/python

"""LM-defined function for strategy description."""

from rdkit.Chem import AllChem, rdFMCS
import copy
from collections import deque
import rdkit.Chem as Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdChemReactions
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
import rdkit.Chem.rdFMCS
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem import AllChem, Descriptors, Lipinski
from rdkit.Chem import rdmolops
import re
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import AllChem, Descriptors
import traceback
import rdkit
from collections import Counter


def main(route):
    """
    Detects if nitrile (C#N) and halide functional groups are preserved
    throughout the synthesis route.
    """
    # Track molecules at each step
    molecules_by_depth = {}
    max_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal max_depth
        max_depth = max(max_depth, depth)

        if node["type"] == "mol" and "smiles" in node:
            if depth not in molecules_by_depth:
                molecules_by_depth[depth] = []
            molecules_by_depth[depth].append(node["smiles"])

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Define SMARTS patterns for functional groups to track
    nitrile_pattern = Chem.MolFromSmarts("[C]#[N]")
    halide_pattern = Chem.MolFromSmarts("[F,Cl,Br,I]")

    # Check if both functional groups are present at each depth
    preserved = True
    for depth in range(max_depth + 1):
        if depth not in molecules_by_depth:
            continue

        depth_has_nitrile = False
        depth_has_halide = False

        for smiles in molecules_by_depth[depth]:
            try:
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    if mol.HasSubstructMatch(nitrile_pattern):
                        depth_has_nitrile = True
                    if mol.HasSubstructMatch(halide_pattern):
                        depth_has_halide = True
            except:
                print(f"Error processing SMILES: {smiles}")

        if not (depth_has_nitrile and depth_has_halide):
            preserved = False
            break

    print(f"Preserved functional groups strategy detected: {preserved}")
    return preserved
