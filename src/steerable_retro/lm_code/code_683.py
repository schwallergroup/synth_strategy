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
    This function detects preservation of stereochemistry throughout the synthesis.
    """
    has_stereocenter = False
    preserves_stereocenter = True

    def dfs_traverse(node):
        nonlocal has_stereocenter, preserves_stereocenter

        if node["type"] == "mol" and "smiles" in node:
            # Check if molecule has stereochemistry
            smiles = node["smiles"]
            if "@" in smiles:
                has_stereocenter = True

                # Count stereocenters
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=False)
                    if len(chiral_centers) == 0:
                        preserves_stereocenter = False
                        print("Found a molecule that lost stereochemistry")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from root
    dfs_traverse(route)

    # Only return true if there was a stereocenter and it was preserved
    if has_stereocenter and preserves_stereocenter:
        print("Stereocenter preservation strategy detected")
        return True
    return False
