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
    Detects if the 1,3-benzodioxole (methylenedioxy) scaffold is preserved throughout the synthesis.
    """
    methylenedioxy_pattern = Chem.MolFromSmarts("c1cc2OCOc2cc1")
    all_steps_have_pattern = True

    def dfs_traverse(node):
        nonlocal all_steps_have_pattern

        if node["type"] == "mol" and node["smiles"]:
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol is not None and not mol.HasSubstructMatch(methylenedioxy_pattern):
                if not node.get("in_stock", False):  # Ignore starting materials
                    all_steps_have_pattern = False
                    print(f"Molecule without methylenedioxy pattern: {node['smiles']}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return all_steps_have_pattern
