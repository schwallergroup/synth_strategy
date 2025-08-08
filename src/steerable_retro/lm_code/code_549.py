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
    Detects if the synthesis route involves carboxamide intermediates
    before final heterocycle formation.
    """
    has_carboxamide_intermediate = False

    def dfs_traverse(node, depth=0):
        nonlocal has_carboxamide_intermediate

        if node["type"] == "mol" and "smiles" in node and depth > 0:
            # Check intermediates (not the final product)
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                carboxamide_pattern = Chem.MolFromSmarts("[#6](=[#8])[#7]")
                if mol.HasSubstructMatch(carboxamide_pattern):
                    has_carboxamide_intermediate = True
                    print(f"Detected carboxamide intermediate at depth {depth}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    return has_carboxamide_intermediate
