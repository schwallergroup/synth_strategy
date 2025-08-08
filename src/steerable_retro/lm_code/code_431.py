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
    This function detects if the synthetic route involves compounds containing
    a trifluoromethyl group.
    """
    contains_trifluoromethyl = False

    def dfs_traverse(node):
        nonlocal contains_trifluoromethyl

        if node["type"] == "mol" and "smiles" in node:
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                trifluoromethyl_pattern = Chem.MolFromSmarts("[c]-[C]([F])([F])[F]")
                if mol.HasSubstructMatch(trifluoromethyl_pattern):
                    print("Found trifluoromethyl group")
                    contains_trifluoromethyl = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return contains_trifluoromethyl
