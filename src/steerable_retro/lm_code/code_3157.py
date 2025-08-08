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
    This function detects a synthetic strategy where an ester functional group
    is preserved throughout the synthesis.
    """
    # Track molecules at each depth
    molecules_by_depth = {}

    # SMARTS pattern for ester
    ester_pattern = Chem.MolFromSmarts("[#6](=[#8])-[#8]-[#6]")

    def dfs_traverse(node, depth=0):
        if node["type"] == "mol":
            # Store molecule at this depth
            if depth not in molecules_by_depth:
                molecules_by_depth[depth] = []

            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                molecules_by_depth[depth].append(mol)

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if ester is present at all depths
    ester_preserved = True
    for depth, mols in molecules_by_depth.items():
        depth_has_ester = False
        for mol in mols:
            if mol.HasSubstructMatch(ester_pattern):
                depth_has_ester = True
                break

        if not depth_has_ester:
            ester_preserved = False
            break

    if ester_preserved:
        print("Detected strategy with ester group preserved throughout synthesis")
    else:
        print("Ester group not preserved throughout the synthesis")

    return ester_preserved
