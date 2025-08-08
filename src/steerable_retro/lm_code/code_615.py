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
    This function detects a strategy involving aniline as a key intermediate
    in multiple steps of the synthesis.
    """
    aniline_count = 0

    def dfs_traverse(node):
        nonlocal aniline_count

        if node["type"] == "mol" and "smiles" in node:
            # Check for aniline pattern
            aniline_pattern = Chem.MolFromSmarts("c-[#7;H1,H2]")
            try:
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol and mol.HasSubstructMatch(aniline_pattern):
                    aniline_count += 1
                    print(f"Detected aniline intermediate: {node['smiles']}")
            except:
                pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # If aniline appears in multiple molecules, it's a key intermediate
    return aniline_count >= 2
