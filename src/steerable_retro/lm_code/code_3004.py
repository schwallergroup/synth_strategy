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
    This function detects if the synthesis utilizes a tosyl group
    as a directing or activating group.
    """
    uses_tosyl = False

    def dfs_traverse(node):
        nonlocal uses_tosyl

        if node["type"] == "mol" and "smiles" in node:
            # Check for tosyl group
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                tosyl_pattern = Chem.MolFromSmarts("c1ccc(S(=O)(=O)[#6])cc1")
                if mol.HasSubstructMatch(tosyl_pattern):
                    print(f"Detected tosyl group in molecule: {node['smiles']}")
                    uses_tosyl = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return uses_tosyl
