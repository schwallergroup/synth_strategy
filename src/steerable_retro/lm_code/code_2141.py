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
    This function detects if a halogen (chlorine) is retained throughout the synthesis.
    """
    halogen_in_final = False
    halogen_in_starting = False

    def dfs_traverse(node, depth=0, is_leaf=False):
        nonlocal halogen_in_final, halogen_in_starting

        if node["type"] == "mol":
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Check for halogen
                halogen_pattern = Chem.MolFromSmarts("[Cl,Br,I,F]")
                has_halogen = mol.HasSubstructMatch(halogen_pattern)

                # If this is the final product (depth 0)
                if depth == 0:
                    halogen_in_final = has_halogen
                    print(f"Halogen in final product: {halogen_in_final}")

                # If this is a starting material (leaf node with no children)
                if node.get("in_stock", False) or (not node.get("children", [])):
                    if has_halogen:
                        halogen_in_starting = True
                        print(f"Halogen in starting material found at depth {depth}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    return halogen_in_final and halogen_in_starting
