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
    This function checks if a nitrile group is maintained throughout the synthesis.
    """
    nitrile_count = 0
    total_reactions = 0

    def dfs_traverse(node):
        nonlocal nitrile_count, total_reactions

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            total_reactions += 1
            rsmi = node["metadata"]["rsmi"]
            product = rsmi.split(">")[-1]

            # Check for nitrile in product
            nitrile_pattern = Chem.MolFromSmarts("[#6]#[#7]")
            try:
                product_mol = Chem.MolFromSmiles(product)
                if product_mol and product_mol.HasSubstructMatch(nitrile_pattern):
                    nitrile_count += 1
                    print("Found nitrile group in reaction product")
            except:
                pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if nitrile is present in most reactions (>80%)
    return total_reactions > 0 and (nitrile_count / total_reactions) > 0.8
