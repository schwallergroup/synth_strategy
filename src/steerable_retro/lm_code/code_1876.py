#!/bin/python

"""LM-defined function for strategy description."""

import copy
import re
from collections import deque

import rdkit
import rdkit.Chem as Chem
from rdkit import Chem
from rdkit.Chem import (
    AllChem,
    Descriptors,
    Lipinski,
    rdChemReactions,
    rdFMCS,
    rdMolDescriptors,
    rdmolops,
)
from rdkit.Chem.Scaffolds import MurckoScaffold


def main(route):
    """
    Detects a linear synthesis strategy where each step involves only one complex fragment
    that is not a starting material.
    """
    is_linear = True
    step_count = 0

    def dfs_traverse(node):
        nonlocal is_linear, step_count

        if node["type"] == "reaction":
            step_count += 1

            # Count complex non-stock fragments
            complex_non_stock_fragments = 0

            for child in node.get("children", []):
                if child["type"] == "mol" and not child.get("in_stock", False):
                    mol = Chem.MolFromSmiles(child["smiles"])
                    if mol and mol.GetNumAtoms() > 10:
                        complex_non_stock_fragments += 1
                        print(f"Found complex non-stock fragment: {child['smiles']}")

            # If more than one complex non-stock fragment, it's not a linear synthesis
            if complex_non_stock_fragments > 1:
                is_linear = False
                print("Found convergent step with multiple complex non-stock fragments")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if it's a linear synthesis with at least 3 steps
    return is_linear and step_count >= 3
