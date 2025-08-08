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
    This function detects if the final product contains both a thiophene ring and an amide bond.
    """
    final_product_has_features = False

    def dfs_traverse(node):
        nonlocal final_product_has_features

        if node["type"] == "mol" and node.get("depth", 0) == 0:  # Final product (depth 0)
            smiles = node.get("smiles", "")
            mol = Chem.MolFromSmiles(smiles)

            if mol:
                # Check for thiophene
                thiophene_pattern = Chem.MolFromSmarts("[#16]1[#6][#6][#6][#6]1")
                # Check for amide
                amide_pattern = Chem.MolFromSmarts("[#6](=O)[#7]")

                if mol.HasSubstructMatch(thiophene_pattern) and mol.HasSubstructMatch(
                    amide_pattern
                ):
                    final_product_has_features = True
                    print("Final product contains both thiophene and amide")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return final_product_has_features
