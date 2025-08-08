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
    This function detects if the synthesis includes N-sulfonylation with a thiophene group.
    """
    has_n_sulfonylation = False

    def dfs_traverse(node):
        nonlocal has_n_sulfonylation

        if node["type"] == "mol" and not node.get("children"):  # Final product
            mol = Chem.MolFromSmiles(node["smiles"])
            if not mol:
                return

            # Pattern for N-SO2-thiophene
            n_sulfonyl_thiophene_pattern = Chem.MolFromSmarts("[#7]S(=O)(=O)c1cccs1")

            if mol.HasSubstructMatch(n_sulfonyl_thiophene_pattern):
                has_n_sulfonylation = True
                print("Found N-sulfonylation with thiophene in final product")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(f"N-sulfonylation with thiophene: {'present' if has_n_sulfonylation else 'absent'}")
    return has_n_sulfonylation
