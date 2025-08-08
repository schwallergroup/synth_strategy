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
    Detects if the synthesis route maintains unprotected phenol groups
    throughout the synthesis without protection/deprotection.
    """
    phenol_at_depths = set()

    def dfs_traverse(node, depth=0):
        if node["type"] == "mol" and "smiles" in node:
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Pattern for phenol
                phenol_pattern = Chem.MolFromSmarts("c[OH]")
                if len(mol.GetSubstructMatches(phenol_pattern)) > 0:
                    phenol_at_depths.add(depth)
                    print(f"Detected phenol at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if phenol is present at multiple depths
    has_unprotected_phenol = len(phenol_at_depths) >= 2

    if has_unprotected_phenol:
        print(f"Confirmed unprotected phenol throughout synthesis at depths: {phenol_at_depths}")

    return has_unprotected_phenol
