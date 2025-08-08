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
    Detects if the synthesis route involves molecules with stereocenters.
    """
    has_stereocenter = False

    def dfs_traverse(node):
        nonlocal has_stereocenter

        if node["type"] == "mol" and "smiles" in node:
            smiles = node["smiles"]
            # Check for @ symbol which indicates stereochemistry in SMILES
            if "@" in smiles:
                has_stereocenter = True
                print("Found molecule with stereocenter")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return has_stereocenter
