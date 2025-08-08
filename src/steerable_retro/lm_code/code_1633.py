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
    This function detects if the synthesis involves a piperazine scaffold.
    """
    has_piperazine_scaffold = False

    def dfs_traverse(node):
        nonlocal has_piperazine_scaffold

        if node["type"] == "mol":
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Check for piperazine pattern
                piperazine_pattern = Chem.MolFromSmarts("[N]1CCN([#6])CC1")
                if mol.HasSubstructMatch(piperazine_pattern):
                    has_piperazine_scaffold = True
                    print("Detected piperazine scaffold in molecule")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return has_piperazine_scaffold
