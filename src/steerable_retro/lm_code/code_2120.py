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
    This function detects if the synthetic route uses TMS-protected alkynes
    as key intermediates.
    """
    has_tms_alkyne = False

    def dfs_traverse(node):
        nonlocal has_tms_alkyne

        if node["type"] == "mol":
            try:
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol:
                    # SMARTS for TMS-alkyne: [C]#[C][Si]([C])([C])[C]
                    tms_alkyne_pattern = Chem.MolFromSmarts("[C]#[C][Si]([C])([C])[C]")
                    if mol.HasSubstructMatch(tms_alkyne_pattern):
                        has_tms_alkyne = True
                        print(f"Detected TMS-protected alkyne: {node['smiles']}")
            except:
                pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return has_tms_alkyne
