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
    Detects if the synthesis route involves Boc-protected nitrogen compounds.
    """
    found_boc_protection = False

    def dfs_traverse(node, depth=0):
        nonlocal found_boc_protection

        # Check both reaction nodes and molecule nodes
        if "smiles" in node:
            smiles = node["smiles"]

            # Check for Boc group
            boc_pattern = Chem.MolFromSmarts("CC(C)(C)OC(=O)N")

            try:
                mol = Chem.MolFromSmiles(smiles)
                if mol and mol.HasSubstructMatch(boc_pattern):
                    print(f"Found Boc-protected compound at depth {depth}")
                    found_boc_protection = True
            except:
                pass

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return found_boc_protection
