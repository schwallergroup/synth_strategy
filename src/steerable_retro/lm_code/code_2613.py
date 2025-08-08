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
    Detects if the synthesis route involves a benzyl ether linker between aromatic rings.
    """
    has_benzyl_ether = False

    def dfs_traverse(node):
        nonlocal has_benzyl_ether

        if node["type"] == "mol":
            # Check for benzyl ether pattern
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                benzyl_ether_pattern = Chem.MolFromSmarts("c-[CH2]-[O]-c")
                if mol.HasSubstructMatch(benzyl_ether_pattern):
                    has_benzyl_ether = True
                    print(f"Detected benzyl ether linker in molecule: {node['smiles']}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return has_benzyl_ether
