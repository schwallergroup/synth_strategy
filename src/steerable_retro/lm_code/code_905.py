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
    This function detects if the route involves the synthesis of a benzochroman scaffold.
    """
    has_benzochroman_scaffold = False

    def dfs_traverse(node):
        nonlocal has_benzochroman_scaffold

        if node["type"] == "mol":
            # Benzochroman scaffold SMARTS pattern
            # This pattern represents the core tricyclic system of benzochroman
            benzochroman_pattern = Chem.MolFromSmarts("[c]1[c][c]2[c]([c][c]1)[O][C][C][C]2")

            mol = Chem.MolFromSmiles(node["smiles"])
            if mol and mol.HasSubstructMatch(benzochroman_pattern):
                has_benzochroman_scaffold = True
                print("Benzochroman scaffold detected")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(f"Benzochroman scaffold synthesis: {has_benzochroman_scaffold}")
    return has_benzochroman_scaffold
