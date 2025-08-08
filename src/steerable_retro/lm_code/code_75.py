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
    This function detects a synthesis that maintains Boc protection throughout
    while building a complex structure.
    """
    boc_maintained = False
    steps_with_boc = 0
    total_steps = 0

    def dfs_traverse(node, depth=0):
        nonlocal steps_with_boc, total_steps

        if node["type"] == "reaction":
            total_steps += 1

            # Extract product
            rsmi = node["metadata"]["rsmi"]
            product_smiles = rsmi.split(">")[-1]
            product_mol = Chem.MolFromSmiles(product_smiles)

            # Check for Boc group
            boc_pattern = Chem.MolFromSmarts("[#6]([#6])([#6])([#6])[#8][#6](=[#8])[#7]")
            if product_mol.HasSubstructMatch(boc_pattern):
                steps_with_boc += 1
                print(f"Detected Boc protection at depth {depth}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Return True if Boc is maintained throughout most of the synthesis
    boc_maintained = (steps_with_boc > 0) and (steps_with_boc / total_steps >= 0.75)
    print(f"Strategy detection result: {boc_maintained}")
    print(f"Steps with Boc: {steps_with_boc}/{total_steps}")

    return boc_maintained
