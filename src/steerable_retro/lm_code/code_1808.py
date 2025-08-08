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
    This function detects if the synthesis incorporates a hydroxymethyl aromatic group.
    """
    hydroxymethyl_aromatic_present = False

    # SMARTS pattern for hydroxymethyl aromatic
    hydroxymethyl_aromatic_pattern = Chem.MolFromSmarts("[c][CH2][OH]")

    def dfs_traverse(node):
        nonlocal hydroxymethyl_aromatic_present

        if node["type"] == "mol":
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol and mol.HasSubstructMatch(hydroxymethyl_aromatic_pattern):
                hydroxymethyl_aromatic_present = True
                print(f"Hydroxymethyl aromatic found in molecule: {node['smiles']}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Hydroxymethyl aromatic incorporation: {hydroxymethyl_aromatic_present}")
    return hydroxymethyl_aromatic_present
