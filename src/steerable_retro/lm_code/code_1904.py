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
    This function detects preservation of key functional groups (cyano, chloro)
    throughout the synthetic route.
    """
    # Track if we find the functional groups in the final product
    cyano_preserved = False
    chloro_preserved = False

    def dfs_traverse(node):
        nonlocal cyano_preserved, chloro_preserved

        if node["type"] == "mol" and node.get("in_stock", False) == False:
            # This is likely the final product
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Check for cyano group
                cyano_pattern = Chem.MolFromSmarts("[C]#[N]")
                if mol.HasSubstructMatch(cyano_pattern):
                    print("Found cyano group in final product")
                    cyano_preserved = True

                # Check for chloro group on aromatic
                chloro_pattern = Chem.MolFromSmarts("[c][Cl]")
                if mol.HasSubstructMatch(chloro_pattern):
                    print("Found chloro group in final product")
                    chloro_preserved = True

        # Check if functional groups are present in starting materials
        if node["type"] == "mol" and node.get("in_stock", False) == True:
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Check for cyano group in starting materials
                cyano_pattern = Chem.MolFromSmarts("[C]#[N]")
                if mol.HasSubstructMatch(cyano_pattern):
                    print("Found cyano group in starting material")

                # Check for chloro group in starting materials
                chloro_pattern = Chem.MolFromSmarts("[c][Cl]")
                if mol.HasSubstructMatch(chloro_pattern):
                    print("Found chloro group in starting material")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if both functional groups are preserved
    return cyano_preserved and chloro_preserved
