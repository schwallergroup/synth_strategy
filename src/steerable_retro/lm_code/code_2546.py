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
    Detects if the synthetic route involves late-stage benzoxazole formation.
    Late stage means in the first half of the synthesis (lower depth numbers).
    """
    benzoxazole_formed = False
    benzoxazole_depth = -1
    max_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal benzoxazole_formed, benzoxazole_depth, max_depth

        if depth > max_depth:
            max_depth = depth

        if node["type"] == "reaction":
            # Extract product SMILES
            rsmi = node["metadata"]["rsmi"]
            product = rsmi.split(">")[-1]

            # Check if product contains benzoxazole
            mol = Chem.MolFromSmiles(product)
            if mol:
                benzoxazole_pattern = Chem.MolFromSmarts("c1nc(c2ccccc2)c(c3ccccc3)o1")
                if mol.HasSubstructMatch(benzoxazole_pattern) and not benzoxazole_formed:
                    benzoxazole_formed = True
                    benzoxazole_depth = depth
                    print(f"Benzoxazole formation detected at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if benzoxazole formation occurs in the first half of the synthesis
    if benzoxazole_formed and benzoxazole_depth <= max_depth / 2:
        print(
            f"Late-stage benzoxazole formation confirmed (depth {benzoxazole_depth} of max {max_depth})"
        )
        return True
    return False
