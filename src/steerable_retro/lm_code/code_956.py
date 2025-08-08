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
    Detects the use of a 3-carbon linker (propanol) to connect two aromatic systems.
    """
    linker_strategy_found = False

    def dfs_traverse(node):
        nonlocal linker_strategy_found

        if node["type"] == "mol":
            try:
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol:
                    # Pattern for 3-carbon linker connecting O and N
                    linker_pattern = Chem.MolFromSmarts("[c]-[#8]-[#6]-[#6]-[#6]-[#7]-c")
                    if mol.HasSubstructMatch(linker_pattern):
                        linker_strategy_found = True
                        print("3-carbon linker strategy detected")
            except:
                pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return linker_strategy_found
