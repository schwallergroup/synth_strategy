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
    Detects if the synthesis route maintains a sulfonyl group throughout the synthesis.
    """
    depths_with_sulfonyl = []

    def dfs_traverse(node, depth=0):
        if "smiles" in node and node["smiles"]:
            # Pattern for sulfonyl group
            sulfonyl_pattern = Chem.MolFromSmarts("[#6][S](=[O])(=[O])[#6]")

            try:
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol and mol.HasSubstructMatch(sulfonyl_pattern):
                    print(f"Sulfonyl group detected at depth {depth}")
                    depths_with_sulfonyl.append(depth)
            except:
                pass

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if sulfonyl group is present at multiple depths
    return len(depths_with_sulfonyl) >= 3  # Present in at least 3 intermediates
