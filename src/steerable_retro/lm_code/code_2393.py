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
    Detects a strategy involving preservation of stereocenters throughout the synthesis.
    """
    stereocenters_by_depth = {}

    def dfs_traverse(node):
        nonlocal stereocenters_by_depth

        if node["type"] == "mol":
            smiles = node["smiles"]
            depth = node.get("depth", 0)

            try:
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    # Find chiral centers
                    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=False)
                    stereocenters_by_depth[depth] = len(chiral_centers)
            except:
                pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Check if stereocenters are preserved (at least 2 stereocenters in final product and maintained throughout)
    if not stereocenters_by_depth:
        return False

    final_stereocenters = stereocenters_by_depth.get(0, 0)
    if final_stereocenters < 2:
        return False

    # Check if the number of stereocenters is maintained or increased throughout
    depths = sorted(stereocenters_by_depth.keys())
    for i in range(len(depths) - 1):
        if stereocenters_by_depth[depths[i]] < stereocenters_by_depth[depths[i + 1]]:
            return False

    print(f"Found stereocenter preservation: {stereocenters_by_depth}")
    return True
