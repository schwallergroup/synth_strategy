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
    This function detects if the synthesis involves both guanidine and sulfonamide
    functional groups in the final product.
    """
    # SMARTS patterns for guanidine and sulfonamide
    guanidine_pattern = Chem.MolFromSmarts("[#6](=[#7])[#7]")
    sulfonamide_pattern = Chem.MolFromSmarts("[#7][#16](=[#8])(=[#8])[#6]")

    contains_both = False

    def dfs_traverse(node):
        nonlocal contains_both

        if node["type"] == "mol" and node.get("smiles"):
            try:
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol:
                    has_guanidine = mol.HasSubstructMatch(guanidine_pattern)
                    has_sulfonamide = mol.HasSubstructMatch(sulfonamide_pattern)

                    if has_guanidine and has_sulfonamide:
                        contains_both = True
                        print("Found both guanidine and sulfonamide groups")
            except:
                print(f"Error processing SMILES: {node['smiles']}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Contains both guanidine and sulfonamide: {contains_both}")
    return contains_both
