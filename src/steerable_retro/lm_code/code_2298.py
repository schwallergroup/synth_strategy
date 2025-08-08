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
    Detects if the synthesis involves CF3-containing compounds throughout.
    """
    # Track if CF3 group is present in the final product
    cf3_in_final = False
    cf3_in_intermediates = 0

    def dfs_traverse(node, depth=0):
        nonlocal cf3_in_final, cf3_in_intermediates

        if node["type"] == "mol":
            # Check for CF3 group
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                cf3_pattern = Chem.MolFromSmarts("[CX4]([FX1])([FX1])[FX1]")
                if mol.HasSubstructMatch(cf3_pattern):
                    if depth == 0:  # Final product
                        cf3_in_final = True
                        print("Found CF3 group in final product")
                    else:
                        cf3_in_intermediates += 1
                        print(f"Found CF3 group in intermediate at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Return True if CF3 is present in final product and at least one intermediate
    return cf3_in_final and cf3_in_intermediates > 0
