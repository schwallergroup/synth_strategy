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
    Detects if the synthesis uses halogenated aromatics (chloro, bromo, iodo)
    as key coupling partners.
    """
    haloarene_count = 0

    def dfs_traverse(node, depth=0):
        nonlocal haloarene_count

        if node["type"] == "mol" and node.get("in_stock", False):
            # Check starting materials for haloarenes
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                haloarene_pattern = Chem.MolFromSmarts("[c][Cl,Br,I]")
                if mol.HasSubstructMatch(haloarene_pattern):
                    haloarene_count += 1
                    print(f"Haloarene detected in starting material: {node['smiles']}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return haloarene_count >= 2  # At least 2 haloarenes used in the synthesis
