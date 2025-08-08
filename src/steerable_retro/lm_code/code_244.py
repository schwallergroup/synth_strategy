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
    This function detects if the synthesis route involves a diarylmethylamine scaffold.
    """
    contains_diarylmethylamine = False

    def dfs_traverse(node):
        nonlocal contains_diarylmethylamine

        if node["type"] == "mol" and "smiles" in node:
            smiles = node["smiles"]

            # Check for diarylmethylamine scaffold
            diarylmethylamine_pattern = Chem.MolFromSmarts(
                "[#7][#6]([#6]1[#6][#6][#6][#6][#6]1)[#6]1[#6][#6][#6][#6][#6]1"
            )
            try:
                mol = Chem.MolFromSmiles(smiles)
                if mol and mol.HasSubstructMatch(diarylmethylamine_pattern):
                    print("Diarylmethylamine scaffold detected")
                    contains_diarylmethylamine = True
            except:
                pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Call dfs_traverse on the root node
    dfs_traverse(route)

    return contains_diarylmethylamine
