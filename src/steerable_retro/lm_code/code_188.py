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
    This function detects if the final product contains multiple heterocycles
    (specifically thiazole, indole, and thiophene).
    """
    has_multiple_heterocycles = False

    def dfs_traverse(node):
        nonlocal has_multiple_heterocycles

        if node["type"] == "mol" and not node.get("children"):  # Final product
            mol = Chem.MolFromSmiles(node["smiles"])
            if not mol:
                return

            # Patterns for heterocycles
            thiazole_pattern = Chem.MolFromSmarts("[#6]1[#16][#6][#6][#7]1")
            indole_pattern = Chem.MolFromSmarts("c1ccc2[nH]ccc2c1")
            thiophene_pattern = Chem.MolFromSmarts("c1cccs1")

            # Count heterocycles
            heterocycle_count = 0

            if mol.HasSubstructMatch(thiazole_pattern):
                heterocycle_count += 1
                print("Found thiazole in final product")

            if mol.HasSubstructMatch(indole_pattern):
                heterocycle_count += 1
                print("Found indole in final product")

            if mol.HasSubstructMatch(thiophene_pattern):
                heterocycle_count += 1
                print("Found thiophene in final product")

            if heterocycle_count >= 2:
                has_multiple_heterocycles = True
                print(f"Found {heterocycle_count} different heterocycles in final product")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(f"Multiple heterocycles: {'present' if has_multiple_heterocycles else 'absent'}")
    return has_multiple_heterocycles
