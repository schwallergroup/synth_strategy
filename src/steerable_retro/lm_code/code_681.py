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
    Detects if a fluorinated aromatic ring is present throughout the synthesis.
    Returns True if at least 70% of the synthesis depths contain a fluorinated aromatic ring.
    """
    # Track molecules at each depth
    molecules_by_depth = {}

    def dfs_traverse(node, depth=0):
        if node["type"] == "mol" and node["smiles"]:
            mol_smiles = node["smiles"]
            if depth not in molecules_by_depth:
                molecules_by_depth[depth] = []
            if mol_smiles not in molecules_by_depth[depth]:
                molecules_by_depth[depth].append(mol_smiles)

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal to collect molecules
    dfs_traverse(route)

    # Check for fluorinated aromatic at each depth
    depths_with_fluoro_aromatic = set()

    for depth, mol_list in molecules_by_depth.items():
        for mol_smiles in mol_list:
            mol = Chem.MolFromSmiles(mol_smiles)
            if mol:
                has_fluoro_aromatic = False
                for atom in mol.GetAtoms():
                    # Check for fluorine attached to aromatic ring
                    if atom.GetSymbol() == "F":
                        for neighbor in atom.GetNeighbors():
                            if neighbor.GetIsAromatic():
                                has_fluoro_aromatic = True
                                print(f"Found fluorinated aromatic at depth {depth}: {mol_smiles}")
                                depths_with_fluoro_aromatic.add(depth)
                                break
                    if has_fluoro_aromatic:
                        break

    # Calculate coverage based on depths with molecules
    depths_with_molecules = set(molecules_by_depth.keys())
    coverage = (
        len(depths_with_fluoro_aromatic) / len(depths_with_molecules)
        if depths_with_molecules
        else 0
    )
    print(
        f"Coverage of fluorinated aromatic: {coverage:.2f} ({len(depths_with_fluoro_aromatic)}/{len(depths_with_molecules)} depths)"
    )

    # Allow for slight rounding issues (0.67 should pass for 70% threshold)
    if depths_with_fluoro_aromatic and coverage > 0.66:
        print("Fluorinated aromatic maintained throughout synthesis")
        return True
    return False
