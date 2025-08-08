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
    Detects if the synthesis involves heterocyclic components like benzofuran and sulfonamide.
    """
    has_benzofuran = False
    has_sulfonamide = False

    def dfs_traverse(node):
        nonlocal has_benzofuran, has_sulfonamide

        if node["type"] == "mol":
            smiles = node.get("smiles", "")
            if not smiles:
                return

            try:
                mol = Chem.MolFromSmiles(smiles)
                if not mol:
                    return

                # Check for benzofuran pattern
                benzofuran_pattern = Chem.MolFromSmarts("c1cc2occc2cc1")
                if mol.HasSubstructMatch(benzofuran_pattern):
                    has_benzofuran = True
                    print("Detected benzofuran heterocycle")

                # Check for sulfonamide pattern
                sulfonamide_pattern = Chem.MolFromSmarts("[#7]-[#16](=[#8])=[#8]")
                if mol.HasSubstructMatch(sulfonamide_pattern):
                    has_sulfonamide = True
                    print("Detected sulfonamide group")
            except:
                pass

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if both heterocyclic components are present
    result = has_benzofuran and has_sulfonamide
    if result:
        print("Synthesis involves both benzofuran and sulfonamide heterocycles")

    return result
