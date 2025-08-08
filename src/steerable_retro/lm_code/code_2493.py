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
    This function detects if a methylthio group is conserved throughout the synthesis.
    """
    methylthio_pattern = Chem.MolFromSmarts("[c][S][CH3]")
    all_mols_have_methylthio = True

    def dfs_traverse(node):
        nonlocal all_mols_have_methylthio

        if node["type"] == "mol" and node.get("in_stock", False) == False:
            # Check if this molecule has a methylthio group
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                if not mol.HasSubstructMatch(methylthio_pattern):
                    all_mols_have_methylthio = False
                    print(f"Molecule without methylthio group found: {node['smiles']}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    print(f"Conserved methylthio group strategy detected: {all_mols_have_methylthio}")
    return all_mols_have_methylthio
