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
    This function detects the use of halogenated aromatics as key intermediates.
    """
    halogenated_aromatics = False

    def dfs_traverse(node):
        nonlocal halogenated_aromatics

        if node["type"] == "mol" and node.get("in_stock", False) == False:
            # Check intermediate molecules for halogenated aromatics
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Check for chlorinated or fluorinated aromatics
                chloro_aromatic = Chem.MolFromSmarts("[c][Cl]")
                fluoro_aromatic = Chem.MolFromSmarts("[c][F]")

                if mol.HasSubstructMatch(chloro_aromatic) or mol.HasSubstructMatch(fluoro_aromatic):
                    print("Halogenated aromatic intermediate detected")
                    halogenated_aromatics = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from root
    dfs_traverse(route)
    return halogenated_aromatics
