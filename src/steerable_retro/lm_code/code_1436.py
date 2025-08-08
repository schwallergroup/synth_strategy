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
    This function detects if the synthetic route involves fluorinated groups,
    particularly trifluoromethyl groups.
    """
    has_fluorinated_groups = False

    def dfs_traverse(node):
        nonlocal has_fluorinated_groups

        if node["type"] == "mol" and "smiles" in node:
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Check for fluorine atoms
                for atom in mol.GetAtoms():
                    if atom.GetAtomicNum() == 9:  # Fluorine
                        has_fluorinated_groups = True
                        break

                # Specifically check for trifluoromethyl group
                trifluoromethyl_pattern = Chem.MolFromSmarts("[CX4]([F])([F])[F]")
                if mol.HasSubstructMatch(trifluoromethyl_pattern):
                    print("Trifluoromethyl group detected")
                    has_fluorinated_groups = True

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from root
    dfs_traverse(route)
    return has_fluorinated_groups
