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
    Detects if the synthesis route involves an enamine intermediate.
    """
    has_enamine = False

    def dfs_traverse(node):
        nonlocal has_enamine

        if node["type"] == "mol":
            # Check for enamine pattern
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                enamine_pattern = Chem.MolFromSmarts("[#6]-[#7](-[#6])-[#6]=[#6]")
                if mol.HasSubstructMatch(enamine_pattern):
                    has_enamine = True
                    print(f"Detected enamine intermediate in molecule: {node['smiles']}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return has_enamine
