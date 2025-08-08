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
    Detects a synthetic strategy where a fluorinated aromatic group is maintained throughout the synthesis.
    """
    fluoro_aromatic_count = 0
    total_reactions = 0

    def dfs_traverse(node):
        nonlocal fluoro_aromatic_count, total_reactions

        if node["type"] == "reaction":
            total_reactions += 1
            rsmi = node["metadata"]["rsmi"]
            product = rsmi.split(">")[-1]

            # Check for fluorinated aromatic in product
            product_mol = Chem.MolFromSmiles(product)
            if product_mol and product_mol.HasSubstructMatch(Chem.MolFromSmarts("c[F]")):
                fluoro_aromatic_count += 1
                print("Found fluorinated aromatic in reaction product")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # Return True if fluorinated aromatic is present in all reactions
    result = fluoro_aromatic_count > 0 and fluoro_aromatic_count == total_reactions
    print(f"Fluorinated aromatic retention strategy detected: {result}")
    return result
