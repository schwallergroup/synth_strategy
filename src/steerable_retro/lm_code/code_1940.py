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
    This function detects if the final product contains a pyridine heterocycle.
    """
    contains_pyridine = False

    def dfs_traverse(node):
        nonlocal contains_pyridine

        if node["type"] == "mol" and not node.get("children"):  # Final product node
            try:
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol:
                    pyridine_pattern = Chem.MolFromSmarts("c1ccncc1")
                    if mol.HasSubstructMatch(pyridine_pattern):
                        print(f"Final product contains pyridine: {node['smiles']}")
                        contains_pyridine = True
            except:
                pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(f"Pyridine in final product: {contains_pyridine}")
    return contains_pyridine
