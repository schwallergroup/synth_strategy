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
    Detects if a TMS-protected alkyne is maintained throughout the synthesis.
    """
    tms_alkyne_depths = []

    def dfs_traverse(node, depth=0):
        nonlocal tms_alkyne_depths

        if node["type"] == "mol":
            smiles = node["smiles"]
            try:
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    # Check for TMS-alkyne
                    tms_alkyne_pattern = Chem.MolFromSmarts("[C]#[C][Si]([C])([C])[C]")
                    if mol.HasSubstructMatch(tms_alkyne_pattern):
                        tms_alkyne_depths.append(depth)
                        print(f"Found TMS-alkyne at depth {depth}: {smiles}")
            except:
                pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if TMS-alkyne appears at multiple depths (maintained throughout synthesis)
    return len(tms_alkyne_depths) >= 2
