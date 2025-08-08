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
    This function detects if the synthesis utilizes isocyanide chemistry
    as a key intermediate step.
    """
    uses_isocyanide = False

    def dfs_traverse(node):
        nonlocal uses_isocyanide

        if node["type"] == "mol" and "smiles" in node:
            # Check for isocyanide group using SMARTS
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                isocyanide_pattern = Chem.MolFromSmarts("[#6]-[#7+]#[#6-]")
                if mol.HasSubstructMatch(isocyanide_pattern):
                    print(f"Detected isocyanide intermediate: {node['smiles']}")
                    uses_isocyanide = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return uses_isocyanide
