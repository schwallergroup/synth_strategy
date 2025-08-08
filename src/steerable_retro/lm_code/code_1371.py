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
    Detects the construction of a pyrrolo[2,3-b]pyridine scaffold.
    """
    scaffold_constructed = False

    def dfs_traverse(node):
        nonlocal scaffold_constructed

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            product = rsmi.split(">")[-1]

            # Check for pyrrolo[2,3-b]pyridine scaffold in product
            scaffold_pattern = Chem.MolFromSmarts("c1cncc2[nH]ccc12")
            try:
                mol = Chem.MolFromSmiles(product)
                if mol and mol.HasSubstructMatch(scaffold_pattern):
                    scaffold_constructed = True
                    print(f"Found pyrrolo[2,3-b]pyridine scaffold in product: {product}")
            except:
                pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return scaffold_constructed
