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
    This function detects a synthetic strategy where a trifluoromethyl group is preserved
    throughout the synthesis.
    """
    trifluoromethyl_present_at_all_steps = True

    def dfs_traverse(node):
        nonlocal trifluoromethyl_present_at_all_steps

        if node["type"] == "mol" and "smiles" in node:
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Check for trifluoromethyl group
                cf3_pattern = Chem.MolFromSmarts("[C]([F])([F])[F]")
                if not mol.HasSubstructMatch(cf3_pattern):
                    # If any molecule doesn't have CF3, set flag to False
                    # But only if it's a significant intermediate, not a reagent
                    if not node.get("in_stock", False) and len(node.get("children", [])) > 0:
                        trifluoromethyl_present_at_all_steps = False
                        print(f"Molecule without CF3 group found: {node['smiles']}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return trifluoromethyl_present_at_all_steps
