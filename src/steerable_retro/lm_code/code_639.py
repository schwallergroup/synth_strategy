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
    Detects synthesis involving a thioether linker connecting heterocyclic systems.
    """
    has_thioether_linker = False
    has_connected_heterocycles = False

    def dfs_traverse(node, depth=0):
        nonlocal has_thioether_linker, has_connected_heterocycles

        if node["type"] == "mol":
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Check for thioether linker
                thioether_linker = Chem.MolFromSmarts("[#6][#16][#6][#6][#7]")
                if mol.HasSubstructMatch(thioether_linker):
                    has_thioether_linker = True
                    print("Detected thioether linker")

                # Check for connected heterocycles
                thiophene = Chem.MolFromSmarts("c1cccs1")
                pyridine = Chem.MolFromSmarts("c1ccncc1")
                pyrimidinone = Chem.MolFromSmarts("c1nc(=O)[nH]cc1")

                if mol.HasSubstructMatch(thiophene) and (
                    mol.HasSubstructMatch(pyridine) or mol.HasSubstructMatch(pyrimidinone)
                ):
                    has_connected_heterocycles = True
                    print("Detected connected heterocycles")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    result = has_thioether_linker and has_connected_heterocycles
    print(f"Heterocycle containing thioether linker: {result}")
    return result
