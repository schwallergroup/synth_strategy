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
    Detects if the synthesis route contains fluoroarene motifs in multiple intermediates.
    """
    fluoroarene_count = 0

    def dfs_traverse(node, depth=0):
        nonlocal fluoroarene_count

        if node["type"] == "mol" and node["smiles"]:
            try:
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol:
                    fluoroarene_pattern = Chem.MolFromSmarts("[c][F]")
                    if mol.HasSubstructMatch(fluoroarene_pattern):
                        print(f"Found fluoroarene motif at depth {depth}")
                        fluoroarene_count += 1
            except:
                pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    return fluoroarene_count >= 3  # Return True if fluoroarene appears in at least 3 intermediates
