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
    This function detects a strategy involving incorporation of a trifluoromethoxy group.
    """
    found_trifluoromethoxy = False

    def dfs_traverse(node, depth=0):
        nonlocal found_trifluoromethoxy

        if node["type"] == "mol":
            # Check if molecule contains trifluoromethoxy group
            smiles = node["smiles"]
            mol = Chem.MolFromSmiles(smiles)

            if mol:
                # SMARTS pattern for trifluoromethoxy group
                trifluoromethoxy_pattern = Chem.MolFromSmarts("[O][C]([F])([F])[F]")

                if mol.HasSubstructMatch(trifluoromethoxy_pattern):
                    found_trifluoromethoxy = True
                    print(f"Found trifluoromethoxy group at depth {depth}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return found_trifluoromethoxy
