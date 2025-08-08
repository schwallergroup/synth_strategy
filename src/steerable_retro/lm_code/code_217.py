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
    This function detects if the final product contains a trifluoromethyl group.
    """
    has_trifluoromethyl = False

    def dfs_traverse(node, depth=0):
        nonlocal has_trifluoromethyl

        if node["type"] == "mol" and depth == 0:  # Final product
            smiles = node["smiles"]
            mol = Chem.MolFromSmiles(smiles)

            # Check for trifluoromethyl group
            trifluoromethyl_pattern = Chem.MolFromSmarts("[#6]([F])([F])[F]")

            if mol and mol.HasSubstructMatch(trifluoromethyl_pattern):
                has_trifluoromethyl = True
                print("Detected trifluoromethyl group in final product")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return has_trifluoromethyl
