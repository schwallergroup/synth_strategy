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
    This function detects a synthesis strategy involving formation of an oxygen-containing
    heterocycle (like tetrahydrofuran) and incorporation of an aryl ether moiety.
    """
    # Track if we found key features
    found_o_heterocycle = False
    found_aryl_ether = False

    # Patterns
    thf_pattern = Chem.MolFromSmarts("[#6]1[#6][#6][#8][#6]1")
    chlorophenoxy_pattern = Chem.MolFromSmarts("[#8][c]1[cH][cH][c]([#17])[cH][cH]1")

    def dfs_traverse(node):
        nonlocal found_o_heterocycle, found_aryl_ether

        if node["type"] == "mol":
            # Check molecule for patterns
            mol = Chem.MolFromSmiles(node["smiles"]) if node["smiles"] else None

            if mol:
                if mol.HasSubstructMatch(thf_pattern):
                    found_o_heterocycle = True
                    print("Oxygen heterocycle (THF) detected")

                if mol.HasSubstructMatch(chlorophenoxy_pattern):
                    found_aryl_ether = True
                    print("Chlorophenoxy group detected")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Check if both features are present
    if found_o_heterocycle and found_aryl_ether:
        print("Oxygen heterocycle with aryl ether strategy detected!")
        return True

    return False
