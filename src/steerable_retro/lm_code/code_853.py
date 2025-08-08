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
    This function detects if a complex heterocyclic core structure is maintained throughout the synthesis.
    """
    # Define a pattern for the heterocyclic core - in this case looking for the benzofuran core
    heterocyclic_core_pattern = Chem.MolFromSmarts("c1ccc2occc2c1")
    core_maintained = True

    def dfs_traverse(node):
        nonlocal core_maintained

        if node["type"] == "mol" and not node.get("in_stock", False):
            try:
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol and not mol.HasSubstructMatch(heterocyclic_core_pattern):
                    core_maintained = False
            except:
                pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    if core_maintained:
        print("Heterocyclic core maintained throughout synthesis")

    return core_maintained
