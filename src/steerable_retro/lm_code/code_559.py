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
    This function detects if the synthesis involves a biphenyl core with
    an ether linker to another functional group.
    """
    biphenyl_pattern = Chem.MolFromSmarts("c1ccccc1-c1ccccc1")
    ether_linker_pattern = Chem.MolFromSmarts("cOCCC")

    contains_biphenyl_with_ether = False

    def dfs_traverse(node):
        nonlocal contains_biphenyl_with_ether

        if node["type"] == "mol" and node["smiles"]:
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                has_biphenyl = mol.HasSubstructMatch(biphenyl_pattern)
                has_ether_linker = mol.HasSubstructMatch(ether_linker_pattern)

                if has_biphenyl and has_ether_linker:
                    contains_biphenyl_with_ether = True
                    print(f"Found biphenyl core with ether linker in molecule: {node['smiles']}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return contains_biphenyl_with_ether
