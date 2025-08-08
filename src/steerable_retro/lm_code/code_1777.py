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
    Detects if the synthetic route maintains a trifluoromethyl group throughout the synthesis.
    """
    has_trifluoromethyl_in_final = False

    def dfs_traverse(node):
        nonlocal has_trifluoromethyl_in_final

        if node["type"] == "mol" and not node.get("children"):
            # This is a leaf node (final product)
            try:
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol:
                    # Check for trifluoromethyl group
                    trifluoromethyl_pattern = Chem.MolFromSmarts("[#6]C(F)(F)F")
                    if mol.HasSubstructMatch(trifluoromethyl_pattern):
                        print("Found trifluoromethyl group in final product")
                        has_trifluoromethyl_in_final = True
            except:
                pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return has_trifluoromethyl_in_final
