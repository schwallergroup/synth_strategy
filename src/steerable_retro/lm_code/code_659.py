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
    This function detects if the synthesis involves a complex heterocyclic scaffold,
    specifically a dibenzothiophene system.
    """
    # SMARTS pattern for dibenzothiophene
    dibenzothiophene_pattern = Chem.MolFromSmarts("c1ccc2c(c1)sc1ccccc12")

    contains_scaffold = False

    def dfs_traverse(node):
        nonlocal contains_scaffold

        if node["type"] == "mol" and node.get("smiles"):
            try:
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol and mol.HasSubstructMatch(dibenzothiophene_pattern):
                    contains_scaffold = True
                    print("Found dibenzothiophene scaffold")
            except:
                print(f"Error processing SMILES: {node['smiles']}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Contains complex heterocyclic scaffold: {contains_scaffold}")
    return contains_scaffold
