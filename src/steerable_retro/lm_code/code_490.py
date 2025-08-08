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
    This function detects a synthetic strategy that incorporates and preserves
    a morpholine-containing aromatic system throughout the synthesis.
    """
    morpholine_pattern = Chem.MolFromSmarts("[#7]1[#6][#6][#8][#6][#6]1")
    morpholine_aromatic_pattern = Chem.MolFromSmarts("[#7]1[#6][#6][#8][#6][#6]1-[c]")

    final_product_has_morpholine = False
    intermediates_with_morpholine = 0

    def dfs_traverse(node, depth=0):
        nonlocal final_product_has_morpholine, intermediates_with_morpholine

        if node["type"] == "mol":
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol and mol.HasSubstructMatch(morpholine_pattern):
                if depth == 0:  # Final product
                    final_product_has_morpholine = True
                    print(f"Final product contains morpholine: {node['smiles']}")
                else:  # Intermediate
                    intermediates_with_morpholine += 1
                    print(f"Intermediate contains morpholine: {node['smiles']}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)

    # Strategy is present if morpholine is in final product and at least one intermediate
    strategy_present = final_product_has_morpholine and intermediates_with_morpholine > 0
    print(f"Morpholine-containing synthesis strategy detected: {strategy_present}")
    return strategy_present
