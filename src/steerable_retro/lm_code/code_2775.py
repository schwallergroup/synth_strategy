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
    This function detects a synthetic strategy involving fluorinated aromatic rings
    in the final product.
    """
    has_fluorinated_aromatic = False

    def dfs_traverse(node, depth=0):
        nonlocal has_fluorinated_aromatic

        if node["type"] == "mol" and depth == 0:  # Check final product
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                fluoro_aromatic_pattern = Chem.MolFromSmarts("[c][F]")
                if mol.HasSubstructMatch(fluoro_aromatic_pattern):
                    print("Found fluorinated aromatic in final product")
                    has_fluorinated_aromatic = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)

    return has_fluorinated_aromatic
