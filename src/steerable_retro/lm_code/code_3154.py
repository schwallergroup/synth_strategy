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
    This function detects a synthetic strategy that involves a diazo compound
    as a key intermediate in the synthesis.
    """
    has_diazo = False

    # SMARTS pattern for diazo group
    diazo_pattern = Chem.MolFromSmarts("[#6]=[#7+]=[#7-]")

    def dfs_traverse(node):
        nonlocal has_diazo

        if node["type"] == "mol":
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol and mol.HasSubstructMatch(diazo_pattern):
                has_diazo = True
                print(f"Found diazo intermediate: {node['smiles']}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    if has_diazo:
        print("Detected synthesis strategy using diazo intermediate")
    else:
        print("No diazo intermediate found in the synthesis route")

    return has_diazo
