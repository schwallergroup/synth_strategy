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
    Detects if the synthesis maintains a stereocenter throughout the route.
    """
    stereocenters_by_depth = {}

    def dfs_traverse(node, depth=0):
        if node["type"] == "mol":
            smiles = node["smiles"]
            try:
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    # Find atoms with specified stereochemistry
                    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=False)
                    if chiral_centers:
                        print(f"Detected {len(chiral_centers)} stereocenters at depth {depth}")
                        stereocenters_by_depth[depth] = len(chiral_centers)
            except:
                print("Error processing molecule SMILES")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if stereocenters are maintained throughout
    depths = sorted(stereocenters_by_depth.keys())
    if len(depths) >= 2 and depths[0] == 0:  # Ensure final product has stereocenter
        print(f"Stereocenters detected at depths: {depths}")
        return True
    return False
