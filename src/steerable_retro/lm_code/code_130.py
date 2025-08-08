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
    final_product_has_cf3 = False

    def dfs_traverse(node, depth=0):
        nonlocal final_product_has_cf3

        if node["type"] == "mol" and depth == 0:  # Final product is at depth 0
            smiles = node["smiles"]
            mol = Chem.MolFromSmiles(smiles)

            if mol is not None:
                cf3_pattern = Chem.MolFromSmarts("FC(F)(F)")
                if mol.HasSubstructMatch(cf3_pattern):
                    print("Final product contains trifluoromethyl group")
                    final_product_has_cf3 = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return final_product_has_cf3
