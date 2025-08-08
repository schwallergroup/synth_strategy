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
    Detects a strategy where a vinyl group is introduced early and preserved throughout the synthesis
    """
    vinyl_depths = []

    def dfs_traverse(node, depth=0):
        nonlocal vinyl_depths

        if node["type"] == "mol":
            # Check for vinyl group
            vinyl_patt = Chem.MolFromSmarts("[#6]=[#6H2]")
            mol = Chem.MolFromSmiles(node["smiles"])

            if mol and mol.HasSubstructMatch(vinyl_patt):
                vinyl_depths.append(depth)
                print(f"Found vinyl group at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if vinyl group appears at multiple depths including high depths (early in synthesis)
    if len(vinyl_depths) >= 2 and max(vinyl_depths) >= 3:
        print("Found vinyl preservation strategy")
        return True

    return False
