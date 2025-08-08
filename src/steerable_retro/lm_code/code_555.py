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
    Detects a synthetic strategy involving halogenated intermediates (bromo or chloro compounds)
    as reactive species.
    """
    has_bromo_intermediate = False
    has_chloro_intermediate = False

    def dfs_traverse(node):
        nonlocal has_bromo_intermediate, has_chloro_intermediate

        if node["type"] == "mol" and "smiles" in node and not node.get("in_stock", False):
            # Only check intermediates, not starting materials
            try:
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol:
                    # Check for halogenated compounds
                    bromo_pattern = Chem.MolFromSmarts("[Br]")
                    chloro_pattern = Chem.MolFromSmarts("[Cl]")

                    if mol.HasSubstructMatch(bromo_pattern):
                        has_bromo_intermediate = True
                        print("Found brominated intermediate")

                    if mol.HasSubstructMatch(chloro_pattern):
                        has_chloro_intermediate = True
                        print("Found chlorinated intermediate")
            except:
                print(f"Error processing molecule SMILES: {node['smiles']}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Strategy is present if either bromo or chloro intermediates are found
    return has_bromo_intermediate or has_chloro_intermediate
