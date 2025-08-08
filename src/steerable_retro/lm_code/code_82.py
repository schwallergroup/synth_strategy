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
    Detects a synthesis route that maintains a methoxymethyl (MOM) protection
    on the indole nitrogen throughout the synthesis.
    """
    has_mom_protection = False

    def dfs_traverse(node, depth=0):
        nonlocal has_mom_protection

        if node["type"] == "mol":
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Check for MOM-protected indole
                mom_indole_pattern = Chem.MolFromSmarts(
                    "[#6]1[#6]2[#7]([CH2]O[CH3])[#6][#6][#6]2[#6][#6][#6]1"
                )
                if mol.HasSubstructMatch(mom_indole_pattern):
                    print("Detected MOM-protected indole")
                    has_mom_protection = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return has_mom_protection
