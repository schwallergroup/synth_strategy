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
    This function detects a synthetic strategy involving trifluoromethyl groups
    in the final product.
    """
    has_trifluoromethyl = False

    def dfs_traverse(node, depth=0):
        nonlocal has_trifluoromethyl

        if node["type"] == "mol" and depth == 0:  # Final product
            smiles = node["smiles"]
            mol = Chem.MolFromSmiles(smiles)

            if mol:
                # Check for trifluoromethyl pattern
                cf3_pattern = Chem.MolFromSmarts("[#6]([#9])([#9])[#9]")
                if mol.HasSubstructMatch(cf3_pattern):
                    match_count = len(mol.GetSubstructMatches(cf3_pattern))
                    print(f"Found {match_count} trifluoromethyl groups in final product")
                    has_trifluoromethyl = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    return has_trifluoromethyl
