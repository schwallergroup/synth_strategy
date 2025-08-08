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
    Detects if the synthesis route involves construction of halogenated aromatic
    heterocycles with both chloro and fluoro substituents.
    """
    has_chloro_aromatic = False
    has_fluoro_aromatic = False
    has_nitrogen_heterocycle = False

    def dfs_traverse(node):
        nonlocal has_chloro_aromatic, has_fluoro_aromatic, has_nitrogen_heterocycle

        if node["type"] == "mol":
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Check for chloro-aromatic
                chloro_aromatic_pattern = Chem.MolFromSmarts("c-[Cl]")
                if mol.HasSubstructMatch(chloro_aromatic_pattern):
                    has_chloro_aromatic = True
                    print("Found chloro-aromatic group")

                # Check for fluoro-aromatic
                fluoro_aromatic_pattern = Chem.MolFromSmarts("c-[F]")
                if mol.HasSubstructMatch(fluoro_aromatic_pattern):
                    has_fluoro_aromatic = True
                    print("Found fluoro-aromatic group")

                # Check for nitrogen heterocycle
                nitrogen_in_ring_pattern = Chem.MolFromSmarts("[#7]@[*]")
                if mol.HasSubstructMatch(nitrogen_in_ring_pattern):
                    has_nitrogen_heterocycle = True
                    print("Found nitrogen heterocycle")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if all conditions are met
    return has_chloro_aromatic and has_fluoro_aromatic and has_nitrogen_heterocycle
