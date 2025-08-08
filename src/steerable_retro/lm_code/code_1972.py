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
    This function detects if the synthesis involves boronic acid/ester chemistry.
    """
    boronic_acid_pattern = Chem.MolFromSmarts("[#6]B([#8])[#8]")
    boronic_ester_pattern = Chem.MolFromSmarts("[#6]B1O[#6][#6]O1")

    boronic_chemistry_used = False

    def dfs_traverse(node):
        nonlocal boronic_chemistry_used

        if node["type"] == "mol":
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                if mol.HasSubstructMatch(boronic_acid_pattern) or mol.HasSubstructMatch(
                    boronic_ester_pattern
                ):
                    boronic_chemistry_used = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(f"Boronic acid chemistry strategy: {boronic_chemistry_used}")
    return boronic_chemistry_used
