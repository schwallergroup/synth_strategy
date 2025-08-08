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
    This function checks if the final product contains multiple chlorinated aromatic rings.
    """
    has_multiple_chloro = False

    def dfs_traverse(node):
        nonlocal has_multiple_chloro

        if node["type"] == "mol" and node.get("in_stock") is False:  # Final product
            smiles = node["smiles"]
            try:
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    # Find all chlorine atoms attached to aromatic carbons
                    chloro_aromatic_pattern = Chem.MolFromSmarts("c[Cl]")
                    matches = mol.GetSubstructMatches(chloro_aromatic_pattern)
                    if len(matches) >= 2:
                        print(f"Found {len(matches)} chlorinated aromatic positions")
                        has_multiple_chloro = True
            except:
                pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return has_multiple_chloro
