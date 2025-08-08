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
    This function detects the use of halogenated aromatic building blocks in the synthesis.
    """
    halogenated_building_blocks = False

    def dfs_traverse(node, depth=0):
        nonlocal halogenated_building_blocks

        if node["type"] == "mol" and node.get("in_stock", False):
            # Check if this starting material has a halogen on an aromatic ring
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol is not None:
                aryl_halide_pattern = Chem.MolFromSmarts("c-[Br,Cl,F,I]")
                if mol.HasSubstructMatch(aryl_halide_pattern):
                    halogenated_building_blocks = True
                    print(f"Found halogenated aromatic building block: {node['smiles']}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return halogenated_building_blocks
