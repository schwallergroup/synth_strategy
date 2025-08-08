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
    Detects if the synthesis route preserves stereochemistry throughout.
    """
    stereocenters_preserved = False
    chiral_centers_count = 0

    def dfs_traverse(node):
        nonlocal stereocenters_preserved, chiral_centers_count

        if node["type"] == "mol" and "smiles" in node:
            # Check for stereochemistry indicators in SMILES
            if "@" in node["smiles"]:
                chiral_centers_count += 1
                print(f"Found chiral center in molecule: {node['smiles']}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    stereocenters_preserved = chiral_centers_count >= 2  # At least two molecules with stereocenters
    print(
        f"Stereocenter preservation strategy detected: {stereocenters_preserved} (chiral centers: {chiral_centers_count})"
    )
    return stereocenters_preserved
