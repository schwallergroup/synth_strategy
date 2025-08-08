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
    This function detects if the synthetic route involves halogen-containing
    intermediates throughout the synthesis.
    """
    has_halogen_intermediates = False
    halogen_count = 0

    def dfs_traverse(node):
        nonlocal has_halogen_intermediates, halogen_count

        if node["type"] == "mol" and "smiles" in node and not node.get("in_stock", False):
            # Check for halogens in intermediates
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                halogen_pattern = Chem.MolFromSmarts("[F,Cl,Br,I]")
                if mol.HasSubstructMatch(halogen_pattern):
                    halogen_count += 1

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # If multiple intermediates contain halogens
    if halogen_count >= 2:
        print(f"Detected {halogen_count} halogen-containing intermediates")
        has_halogen_intermediates = True

    return has_halogen_intermediates
