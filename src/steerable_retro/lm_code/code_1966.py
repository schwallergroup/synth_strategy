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
    This function detects if a carboxylic acid group is preserved throughout the synthesis.
    """
    steps_with_carboxylic_acid = 0
    total_steps = 0

    def dfs_traverse(node, depth=0):
        nonlocal steps_with_carboxylic_acid, total_steps

        if node["type"] == "reaction":
            total_steps += 1
            rsmi = node["metadata"]["rsmi"]
            product = rsmi.split(">")[-1]

            # Check if product contains carboxylic acid
            carboxylic_acid_pattern = Chem.MolFromSmarts("[#6]C(=O)[OH]")
            prod_mol = Chem.MolFromSmiles(product)

            if prod_mol and prod_mol.HasSubstructMatch(carboxylic_acid_pattern):
                steps_with_carboxylic_acid += 1
                print(f"Carboxylic acid found at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # If carboxylic acid is present in all steps, it's preserved
    return total_steps > 0 and steps_with_carboxylic_acid == total_steps
