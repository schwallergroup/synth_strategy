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
    This function detects if a tetrazole moiety is present throughout the entire synthesis route.
    """
    tetrazole_pattern = Chem.MolFromSmarts("[#7]1[#7][#7][#7]C1=O")
    tetrazole_present_at_all_depths = {}

    def dfs_traverse(node, depth=0):
        if node["type"] == "mol":
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol and mol.HasSubstructMatch(tetrazole_pattern):
                tetrazole_present_at_all_depths[depth] = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if tetrazole is present at all depths
    max_depth = (
        max(tetrazole_present_at_all_depths.keys()) if tetrazole_present_at_all_depths else -1
    )
    all_depths_present = all(
        depth in tetrazole_present_at_all_depths for depth in range(max_depth + 1)
    )

    print(f"Tetrazole preservation check: {all_depths_present}")
    return all_depths_present
