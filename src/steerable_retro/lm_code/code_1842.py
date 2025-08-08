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
    This function detects if the synthesis involves compounds with methoxy groups on aromatic rings.
    """
    methoxyphenyl_count = 0

    def dfs_traverse(node):
        nonlocal methoxyphenyl_count

        if node["type"] == "mol" and "smiles" in node:
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Count methoxy groups on aromatic rings
                matches = mol.GetSubstructMatches(Chem.MolFromSmarts("[c][O][C]"))
                if matches:
                    methoxyphenyl_count = max(methoxyphenyl_count, len(matches))
                    print(f"Found compound with {len(matches)} methoxy groups")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if we found multiple (≥2) methoxy groups
    return methoxyphenyl_count >= 2
