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
    This function detects whether stereochemistry is preserved throughout the synthesis.
    """
    # Track molecules with stereocenters at each depth
    stereocenters_by_depth = {}

    def dfs_traverse(node, depth=0):
        if node["type"] == "mol":
            try:
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol:
                    # Check for stereocenters
                    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=False)
                    if chiral_centers:
                        print(f"Stereocenter(s) found at depth {depth}: {node['smiles']}")
                        stereocenters_by_depth[depth] = chiral_centers
            except:
                print(f"Error processing molecule SMILES: {node['smiles']}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Strategy is present if stereocenters are maintained across multiple depths
    strategy_present = len(stereocenters_by_depth) >= 2

    print(f"Stereochemistry preservation strategy detected: {strategy_present}")
    return strategy_present
