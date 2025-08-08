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
    This function detects if the synthesis involves building around a piperazine scaffold.
    """
    has_piperazine = False
    has_substituted_piperazine = False

    def dfs_traverse(node, depth=0):
        nonlocal has_piperazine, has_substituted_piperazine

        if node["type"] == "mol":
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Basic piperazine pattern
                piperazine_pattern = Chem.MolFromSmarts("[N]1CCN([H,C,c])CC1")

                # Disubstituted piperazine pattern
                disubst_piperazine_pattern = Chem.MolFromSmarts("[N]1CCN([C,c])CC1")

                if mol.HasSubstructMatch(piperazine_pattern):
                    has_piperazine = True

                if mol.HasSubstructMatch(disubst_piperazine_pattern):
                    has_substituted_piperazine = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Return True if we find both basic piperazine and substituted piperazine
    # indicating the synthesis involves building around this scaffold
    return has_piperazine and has_substituted_piperazine
