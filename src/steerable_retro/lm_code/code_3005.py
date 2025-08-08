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
    This function detects the use of multiple protecting groups (phthalimide and Boc)
    in the same synthetic route.
    """
    # Track presence of protecting groups
    phthalimide_present = False
    boc_present = False

    # Define SMARTS patterns
    phthalimide_pattern = Chem.MolFromSmarts("[#7]1C(=O)c2ccccc2C1=O")
    boc_pattern = Chem.MolFromSmarts("[#7]C(=O)OC(C)(C)C")

    def dfs_traverse(node):
        nonlocal phthalimide_present, boc_present

        if node["type"] == "mol":
            try:
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol:
                    # Check for phthalimide group
                    if mol.GetSubstructMatches(phthalimide_pattern):
                        phthalimide_present = True
                        print(f"Phthalimide group found in molecule: {node['smiles']}")

                    # Check for Boc group
                    if mol.GetSubstructMatches(boc_pattern):
                        boc_present = True
                        print(f"Boc group found in molecule: {node['smiles']}")
            except:
                print(f"Error processing molecule SMILES: {node['smiles']}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Strategy is present if both protecting groups are used
    strategy_present = phthalimide_present and boc_present

    print(f"Multiple protecting groups strategy detected: {strategy_present}")
    return strategy_present
