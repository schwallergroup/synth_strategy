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
    This function detects if a nitrile group is maintained throughout the synthesis.
    It checks if the final product and at least one early intermediate contain a nitrile group.
    """
    # Track if we've found nitrile in final product and early intermediate
    final_product_has_nitrile = False
    early_intermediate_has_nitrile = False
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal final_product_has_nitrile, early_intermediate_has_nitrile, max_depth

        # Update max depth
        max_depth = max(max_depth, depth)

        if node["type"] == "mol":
            nitrile_pattern = Chem.MolFromSmarts("[C]#[N]")
            mol = Chem.MolFromSmiles(node["smiles"]) if node["smiles"] else None

            if mol and mol.HasSubstructMatch(nitrile_pattern):
                if depth == 0:  # Final product
                    final_product_has_nitrile = True
                    print("Final product has nitrile group")
                elif depth >= 3:  # Early intermediate (depth >= 3)
                    early_intermediate_has_nitrile = True
                    print(f"Early intermediate at depth {depth} has nitrile group")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Return True if nitrile is maintained from early intermediate to final product
    return final_product_has_nitrile and early_intermediate_has_nitrile and max_depth >= 3
