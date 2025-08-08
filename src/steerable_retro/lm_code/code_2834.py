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
    This function detects a synthetic strategy involving a hydroxylated piperidine intermediate.
    """
    # Track if we found the hydroxylated piperidine
    found_hydroxylated_piperidine = False

    def dfs_traverse(node):
        nonlocal found_hydroxylated_piperidine

        if node["type"] == "mol" and node.get("smiles"):
            # Check for hydroxylated piperidine
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                hydroxylated_piperidine_pattern = Chem.MolFromSmarts(
                    "N1[CH2][CH2][CH]([OH])[CH2][CH2]1"
                )
                if mol.HasSubstructMatch(hydroxylated_piperidine_pattern):
                    print(f"Found hydroxylated piperidine intermediate: {node['smiles']}")
                    found_hydroxylated_piperidine = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return found_hydroxylated_piperidine
