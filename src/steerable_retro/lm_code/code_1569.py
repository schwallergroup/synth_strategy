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
    Detects a synthetic strategy where a nitrile group is present
    throughout the synthesis and retained in the final product.
    """
    nitrile_in_final_product = False
    nitrile_in_intermediates = False

    def dfs_traverse(node, depth=0):
        nonlocal nitrile_in_final_product, nitrile_in_intermediates

        if node["type"] == "mol" and "smiles" in node:
            smiles = node["smiles"]
            mol = Chem.MolFromSmiles(smiles)

            if mol and mol.HasSubstructMatch(Chem.MolFromSmarts("C#N")):
                if depth == 0:  # Final product
                    nitrile_in_final_product = True
                    print("Found nitrile in final product")
                else:  # Intermediate
                    nitrile_in_intermediates = True
                    print("Found nitrile in intermediate")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if nitrile is retained throughout synthesis
    return nitrile_in_final_product and nitrile_in_intermediates
