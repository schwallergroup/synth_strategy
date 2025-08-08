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
    This function detects whether a thiazole-pyrimidine heterocyclic scaffold
    is preserved throughout the synthesis.
    """
    scaffold_preserved = True
    thiazole_pattern = Chem.MolFromSmarts("[#6]1[#6][#16][#6][#6]1")
    pyrimidine_pattern = Chem.MolFromSmarts("[#6]1[#7][#6][#7][#6][#6]1")

    def dfs_traverse(node):
        nonlocal scaffold_preserved

        if node["type"] == "mol" and "smiles" in node:
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                has_thiazole = mol.HasSubstructMatch(thiazole_pattern)
                has_pyrimidine = mol.HasSubstructMatch(pyrimidine_pattern)

                # If this is a significant molecule (not just a reagent)
                # and it doesn't have both scaffolds, mark as not preserved
                if not (has_thiazole and has_pyrimidine) and not node.get("in_stock", False):
                    scaffold_preserved = False
                    print(f"Scaffold not preserved in molecule: {node['smiles']}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    if scaffold_preserved:
        print("Heterocyclic scaffold preserved throughout synthesis")

    return scaffold_preserved
