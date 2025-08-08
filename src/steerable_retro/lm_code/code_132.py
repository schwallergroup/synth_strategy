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
    Detects if the synthesis uses an azide intermediate in the route.
    """
    # Initialize tracking variable
    has_azide_intermediate = False

    # SMARTS pattern for azide
    azide_pattern = Chem.MolFromSmarts("[N-]=[N+]=[N]")

    def dfs_traverse(node):
        nonlocal has_azide_intermediate

        if node["type"] == "mol" and not node.get("in_stock", False):
            # Check if this molecule is an intermediate (not a starting material)
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol and mol.HasSubstructMatch(azide_pattern):
                has_azide_intermediate = True
                print("Found azide intermediate:", node["smiles"])

        for child in node.get("children", []):
            dfs_traverse(child)

    # Traverse the route
    dfs_traverse(route)

    return has_azide_intermediate
