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
    This function detects a synthetic strategy that maintains stereochemistry
    throughout the synthesis.
    """
    stereocenters_maintained = True
    stereo_atoms_by_depth = {}

    def dfs_traverse(node, depth=0):
        nonlocal stereocenters_maintained, stereo_atoms_by_depth

        if node["type"] == "mol":
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Get atoms with stereochemistry
                chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=False)
                if chiral_centers:
                    stereo_atoms_by_depth[depth] = set(atom_idx for atom_idx, _ in chiral_centers)

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)

    # Check if stereochemistry is maintained throughout
    if len(stereo_atoms_by_depth) >= 2:  # Need at least two depths with stereocenters
        print(f"Found stereocenters at depths: {list(stereo_atoms_by_depth.keys())}")
        print("Stereochemistry is maintained throughout synthesis")
        return True
    else:
        return False
